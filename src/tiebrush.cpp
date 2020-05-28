#include "GSam.h"
#include "GArgs.h"
#include "tmerge.h"
#include "GBitVec.h"

#include <set>
#include <functional>
#include <iostream>
#include <map>

#define VERSION "0.0.1"

const char* USAGE="TieBrush v" VERSION " usage:\n"
" tiebrush [-o <outfile>.bam] list.txt | in1.bam in2.bam ...\n"
" Other options: \n"
"  --cigar,-C   : merge if only CIGAR string is the same\n"
"  --clip,-P    : merge if clipped CIGAR string is the same\n"
"  --exon,-E    : merge if exon boundaries are the same\n"
"  --bedcov,-B  : report coverage information in BED format\n";

enum TMrgStrategy {
	tMrgStratFull=0, // same CIGAR and MD
	tMrgStratCIGAR,  // same CIGAR (MD may differ)
	tMrgStratClip,   // same CIGAR after clipping
	tMrgStratExon    // same exons
};

TMrgStrategy mrgStrategy=tMrgStratFull;
TInputFiles inRecords;

GStr outfname;
GSamWriter* outfile=NULL;

uint64_t inCounter=0;
uint64_t outCounter=0;

bool debugMode=false;
bool verbose=false;

bool report_bedcov=false;

int cmpFull(GSamRecord& a, GSamRecord& b) {
	//-- CIGAR && MD strings
	if (a.get_b()->core.n_cigar!=b.get_b()->core.n_cigar) return ((int)a.get_b()->core.n_cigar - (int)b.get_b()->core.n_cigar);
	int cigar_cmp=0;
	if (a.get_b()->core.n_cigar>0)
		cigar_cmp=memcmp(bam_get_cigar(a.get_b()) , bam_get_cigar(b.get_b()), a.get_b()->core.n_cigar*sizeof(uint32_t) );
	if (cigar_cmp!=0) return cigar_cmp;
	// compare MD tag
	char* aMD=a.tag_str("MD");
	char* bMD=b.tag_str("MD");
	if (aMD==NULL || bMD==NULL) {
		if (aMD==bMD) return 0;
		if (aMD!=NULL) return 1;
		return -1;
	}
    return strcmp(aMD, bMD);
}

int cmpCigar(GSamRecord& a, GSamRecord& b) {
	if (a.get_b()->core.n_cigar!=b.get_b()->core.n_cigar) return ((int)a.get_b()->core.n_cigar - (int)b.get_b()->core.n_cigar);
	if (a.get_b()->core.n_cigar==0) return 0;
	return memcmp(bam_get_cigar(a.get_b()) , bam_get_cigar(b.get_b()), a.get_b()->core.n_cigar*sizeof(uint32_t) );
}

int cmpCigarClip(GSamRecord& a, GSamRecord& b) {
	uint32_t a_clen=a.get_b()->core.n_cigar;
	uint32_t b_clen=b.get_b()->core.n_cigar;
	uint32_t* a_cstart=bam_get_cigar(a.get_b());
	uint32_t* b_cstart=bam_get_cigar(b.get_b());
	while (a_clen>0 &&
			((*a_cstart) & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { a_cstart++; a_clen--; }
	while (a_clen>0 &&
			(a_cstart[a_clen-1] & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) a_clen--;
	while (b_clen>0 &&
			((*b_cstart) & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { b_cstart++; b_clen--; }
	while (b_clen>0 &&
			(b_cstart[b_clen-1] & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) b_clen--;
	if (a_clen!=b_clen) return (int)a_clen-(int)b_clen;
	if (a_clen==0) return 0;
	return memcmp(a_cstart, b_cstart, a_clen*sizeof(uint32_t));
}

int cmpExons(GSamRecord& a, GSamRecord& b) {
	if (a.exons.Count()!=b.exons.Count()) return (a.exons.Count()-b.exons.Count());
	for (int i=0;i<a.exons.Count();i++) {
		if (a.exons[i].start!=b.exons[i].start)
			return ((int)a.exons[i].start-(int)b.exons[i].start);
		if (a.exons[i].end!=b.exons[i].end)
			return ((int)a.exons[i].end-(int)b.exons[i].end);
	}
	return 0;
}


//keep track of all SAM alignments starting at the same coordinate
// that were merged into a single alignment
class SPData {
    bool settled; //real SP data, owns its r data and deallocates it on destroy
  public:
	int64_t accYC; //if any merged records had YC tag, their values are accumulated here
	int64_t accYX; //if any merged records had YX tag, their values are accumulated here
	GBitVec* samples; //which samples were the collapsed ones coming from
	                 //number of bits set will be stored as YX:i:(samples.count()+accYX)
	int dupCount; //duplicity count - how many single-alignments were merged into r
	              // will be stored as tag YC:i:(dupCount+accYC)
	GSamRecord* r;
	char tstrand; //'-','+' or '.'
    SPData(GSamRecord* rec=NULL):settled(false), accYC(0), accYX(0), samples(NULL),
    		dupCount(0), r(rec), tstrand('.') {
    	if (r!=NULL) tstrand=r->spliceStrand();
    }

    ~SPData() {
    	if (settled && r!=NULL) delete r;
    	if (samples!=NULL) delete samples;
    }

    void detach(bool dontFree=true) { settled=!dontFree; }
    //detach(true) must never be called before settle()

    void settle(int sampleIdx=-1) { //becomes a standalone SPData record
    	// duplicates the current record
    	if (!settled) { //r is not owned
    		settled=true;
    		GSamRecord* rdup=new GSamRecord(*r);
    		r=rdup;
    	}
    	if (samples==NULL) samples=new GBitVec(inRecords.readers.Count());
    	accYC=r->tag_int("YC");
    	if (accYC==0) ++dupCount;
    	accYX=r->tag_int("YX");
    	if (accYX==0 && sampleIdx>=0) samples->set(sampleIdx);
    }

    void dupAdd(GSamRecord* rec, int sampleIdx=-1) { //merge an external SAM record into this one
    	if (!settled) GError("Error: cannot merge a duplicate into a non-settled SP record!\n");
    	//WARNING: rec MUST be a "duplicate" of current record r
    	int64_t rYC=rec->tag_int("YC");
    	if (rYC) accYC+=rYC;
    	   else ++dupCount;
    	int64_t rYX=rec->tag_int("YX");
    	if (rYX) accYX+=rYX;
    	   else if (sampleIdx>=0) samples->set(sampleIdx);
    }

    bool operator<(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId()) return (r->refId()<b.r->refId());
    	//NOTE: already assuming that start&end must match, no matter the merge strategy
    	if (r->start!=b.r->start) return (r->start<b.r->start);
    	if (tstrand!=b.tstrand) return (tstrand<b.tstrand);
    	if (r->end!=b.r->end) return (r->end<b.r->end);
    	if (mrgStrategy==tMrgStratFull) {
    		int ret=cmpFull(*r, *(b.r));
    		return (ret<0);
    	}
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b.r))<0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b.r))<0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b.r))<0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }

    bool operator==(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId() || r->start!=b.r->start || tstrand!=b.tstrand ||
    			r->end!=b.r->end) return false;
    	if (mrgStrategy==tMrgStratFull) return (cmpFull(*r, *(b.r))==0);
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b.r))==0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b.r))==0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b.r))==0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }
};

void processOptions(int argc, char* argv[]);

void addPData(TInputRecord& irec, GList<SPData>& spdlst) {
  //add and collapse if match found
	SPData* newspd=new SPData(irec.brec);
	if (spdlst.Count()>0) {
		//find if irec can merge into existing SPData
		SPData* spf=spdlst.AddIfNew(newspd, false);
		if (spf!=newspd) { //matches existing SP entry spf
			spf->dupAdd(irec.brec, irec.fidx); //update existing SP entry
			delete newspd;
			return;
		} // not a novel SP data
		//else newspd was added as a separate entry of spdlst
	}
	else { // empty list, just add this
		spdlst.Add(newspd);
	}
	newspd->settle(irec.fidx); //keep its own SAM record copy
}

std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> > cov_pq; // since RB+ Tree is sorted by std::greater - smallest position (pair of seqid and pos) will always be at the top (begin())
std::pair<std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> >::iterator,bool>  cpq_it;

struct CovInterval{
    int tid=-1;
    int start=-1; // start position
    int end=-1; // end position
    int cov=-1; // current coverage for the interval

    void extend(bam_hdr_t *al_hdr,int new_tid,int pos,int new_cov){
        if(tid == -1){
            tid=new_tid;
            start=pos;
            end=pos;
            cov=new_cov;
        }
        else if(tid==new_tid && pos-end==1 && cov==new_cov){ // position can extend the interval
            end++;
        }
        else{
            write(al_hdr);
            tid=new_tid;
            start=pos;
            end=pos;
            cov=new_cov;
        }
    }

    void write(sam_hdr_t *al_hdr){
        fprintf(stdout,"%s\t%d\t%d\t%d\n",al_hdr->target_name[tid],start,end,cov);
    }
} covInterval;

void cleanPriorityQueue(sam_hdr_t* al_hdr){ // cleans all entries
    for(auto& pq : cov_pq){
        covInterval.extend(al_hdr,pq.first.first,pq.first.second,pq.second);
    }
    cov_pq.clear();
}

void cleanPriorityQueue(sam_hdr_t* al_hdr,int tid,int pos){ // process and remove all positions upto but not including the given position
    if(cov_pq.empty()){
        return;
    }
    while(!cov_pq.empty() && cov_pq.begin()->first < std::pair<uint16_t,uint32_t>(tid,pos)){
        covInterval.extend(al_hdr,cov_pq.begin()->first.first,cov_pq.begin()->first.second,cov_pq.begin()->second);
        cov_pq.erase(cov_pq.begin());
    }
}

void incCov(sam_hdr_t *al_hdr,bam1_t *in_rec,int dupCount){
    // remove all elements from the priority queue which indicate positions smaller than current
    // iterate over positions following cigar and add to the priority queue
    int pos=in_rec->core.pos+1; // 1-based

    bool first_match = true; // indicates whether the first match position has been encountered - triggers cleanup of the priority queue

    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);

        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CMATCH: // skip to the next position - no change in coverage
                if(first_match){
                    cleanPriorityQueue(al_hdr,in_rec->core.tid,pos);
                    first_match=false;
                }
                for(int i=0;i<oplen;i++){ // add coverage for each position
                    cpq_it = cov_pq.insert(std::make_pair(std::make_pair(in_rec->core.tid,pos),0));
                    cpq_it.first->second+=dupCount;
                    pos+=1;
                }
                break;
            default:
                fprintf(stderr,"ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
                exit(-1);
        }
    }
}

void flushPData(GList<SPData>& spdlst){ //write spdata to outfile
  if (spdlst.Count()==0) return;
  // write SAM records in spdata to outfile
  for (int i=0;i<spdlst.Count();++i) {
	  SPData& spd=*(spdlst.Get(i));
	  int64_t accYC=spd.accYC+spd.dupCount;
	  int64_t accYX=spd.accYX;
	  if (spd.dupCount>1) { //just to save an unnecessary GBitVec::count() call?
	   accYX+=spd.samples->count();
	  }
	  if (accYC>1) spd.r->add_int_tag("YC", accYC);
	  if (accYX>1) spd.r->add_int_tag("YX", accYX);
      if(report_bedcov){incCov(outfile->get_header(),spd.r->get_b(),accYC);}
	  outfile->write(spd.r);
	  outCounter++;
  }

  spdlst.Clear();
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	inRecords.start();
	if (outfname.is_empty()) outfname="-";
	GSamFileType oftype=(outfname=="-") ?
			GSamFile_SAM : GSamFile_BAM;
	outfile=new GSamWriter(outfname, inRecords.header(), oftype);

	 TInputRecord* irec=NULL;
	 GSamRecord* brec=NULL;
	 GList<SPData> spdata(true, true, true); //list of Same Position data, with all possibly merged records
	 //bool newChr=false;
	 int prev_pos=-1;
	 int prev_tid=-1;
	 while ((irec=inRecords.next())!=NULL) {
		 brec=irec->brec;
		 if (brec->isUnmapped()) continue;
		 inCounter++;
		 int tid=brec->refId();
		 int pos=brec->start; //1-based
		 if (tid!=prev_tid) {
		     if(report_bedcov){cleanPriorityQueue(outfile->get_header());};
			 prev_tid=tid;
			 //newChr=true;
			 prev_pos=-1;
		 }
		 if (pos!=prev_pos) { //new position
		     // all possible reads collapsed for the current set of coordinates
		     // add each position to the priority queue. If already exists - add the current count of reads in the collapsed record
		     flushPData(spdata);
			 prev_pos=pos;
		 }
		 addPData(*irec, spdata);
	 }
    flushPData(spdata);
    if(report_bedcov){cleanPriorityQueue(outfile->get_header());}
    if(report_bedcov){covInterval.write(outfile->get_header());}
	delete outfile;
	inRecords.stop();
    //if (verbose) {
    	double p=100.00 - (double)(outCounter*100.00)/(double)inCounter;
    	GMessage("%ld input records written as %ld (%.2f%% reduction)\n", inCounter, outCounter, p);
    //}
}
// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
	GArgs args(argc, argv, "help;debug;verbose;version;cigar;clip;exon;bedcov;CPEBDVho:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
		GMessage(USAGE);
		exit(1);
	}

	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	bool stratC=(args.getOpt("cigar")!=NULL || args.getOpt('C')!=NULL);
	bool stratP=(args.getOpt("clip")!=NULL || args.getOpt('P')!=NULL);
	bool stratE=(args.getOpt("exon")!=NULL || args.getOpt('E')!=NULL);
    report_bedcov=(args.getOpt("bedcov")!=NULL || args.getOpt('B')!=NULL);
	if (stratC | stratP | stratE) {
		if (!(stratC ^ stratP ^ stratE))
			GError("Error: only one merging strategy can be requested.\n");
		if (stratC) mrgStrategy=tMrgStratCIGAR;
		else if (stratP) mrgStrategy=tMrgStratClip;
		else mrgStrategy=tMrgStratExon;
	}

	debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n", VERSION);
	   exit(0);
	}
	 //verbose=(args.getOpt('v')!=NULL);
	if (verbose) {
	   fprintf(stderr, "Running TieBrush " VERSION ". Command line:\n");
	   args.printCmdLine(stderr);
	}
	outfname=args.getOpt('o');
	const char* ifn=NULL;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		inRecords.Add(ifn);
	}
}


