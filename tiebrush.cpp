#include "GSam.h"
#include "GArgs.h"
#include "tmerge.h"
#include "GBitVec.h"

#define VERSION "0.0.2"

const char* USAGE="TieBrush v" VERSION " usage:\n"
" tiebrush [-o <outfile>.bam] list.txt | in1.bam in2.bam ...\n"
" Other options: \n"
"  --cigar,-C   : merge if only CIGAR string is the same\n"
"  --clip,-P    : merge if clipped CIGAR string is the same\n"
"  --exon,-E   : merge if exon boundaries are the same\n";

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

struct GSegNode {
	GSeg seg;
	GSegNode* next;
	inline uint start() { return seg.start; }
	inline uint end() { return seg.end; }
	GSegNode(uint segstart=0, uint segend=0, GSegNode* nnext=NULL):seg(segstart, segend),
			next(nnext) {}
	GSegNode(GSeg& exon, GSegNode* nnext=NULL):seg(exon),
			next(nnext) {}
};

struct GSegList { //per sample per strand
  GSegNode* startNode;
  uint last_pos;
  int last_dist;
  GSegList():startNode(NULL),last_pos(0),last_dist(-1) {
	  //a new list always has (0,0) interval -- no longer needed
	  //startNode=new GSegNode();
  }

  ~GSegList() {
	  clear();
  }

  void reset() {
	  clear();
	  last_pos=0;
	  last_dist=-1;
	  //startNode=new GSegNode();
	  startNode=NULL;
  }

  void clear() { //delete all nodes
	  GSegNode* p=startNode;
	  GSegNode* next;
	  while (p) {
		  next=p->next;
		  delete p;
		  p=next;
	  }
	  startNode=NULL;
  }

  void clearTo(GSegNode* toNode) {
	  //clear every node up to and *including* toNode
	  GSegNode* p=startNode;
	  GSegNode* next;
	  while (p && p!=toNode) {
		  next=p->next;
		  delete p;
		  p=next;
	  }
	  if (p==NULL) GError("Error: clearTo did not find target node %d-%d!\n",
			  toNode->start(), toNode->end());
	  next=toNode->next;
	  delete toNode;
	  startNode=next;
  }

 void mergeRead(GSamRecord& r) {
	 if (startNode==NULL) {
		 startNode=new GSegNode(r.exons[0]);
		 GSegNode* cn=startNode;
		 for (int i=1;i<r.exons.Count();i++){
			 GSegNode* n=new GSegNode(r.exons[i]);
			 cn->next=n;
			 cn=n;
		 }
		 return;
	 }
	 GSegNode *n=startNode;
	 GSegNode *prev=NULL;
	 for (int i=0;i<r.exons.Count();i++) {
		 GSeg& e=r.exons[i];
		 while (n) {
            if (e.end < n->start()) {
              //exon should be inserted before n!
              GSegNode* nw=new GSegNode(e, n);
              if (n==startNode)
            	startNode=nw;
              else
                prev->next=nw;
              //n is unchanged
          	  prev=nw;
              break; //check next exon against n
            }
            // e.end >= node start
            if (e.start <= n->end()) { //overlap!
              //union should replace/change n
              if (e.start<n->start()) n->seg.start=e.start;
              if (e.end>n->end()) n->seg.end=e.end;
              //we have to check if next nodes are now overlapped as well
              GSegNode* next=n->next;
              while (next && next->start()<=n->end()) {
            	  //overlap, swallow this next node
            	  uint nend=next->end();
            	  n->next=next->next;
            	  delete next;
            	  next=n->next;
                  if (nend>n->end()) {
                	n->seg.end=nend;
                	break; //there cannot be overlap with next node
                  }
              } //while overlap with next nodes
              break; //check next exon against this extended n
            }
            // e.start > node end
            prev=n;
            n=n->next;
		 }
	 }
 }

 int processRead(GSamRecord& r) { //return d=current bundle extent upstream
	 //if the read starts after a gap, d=0
	 //this should only be called ONCE per collapsed read and sample
	 //should NOT be called on reads coming from already merged samples!
	 if (last_pos==r.start) { //already called on the same sample and start position
	     mergeRead(r);
		 return last_dist;
	 }
	 int d=0;
	 GSegNode* node=startNode;
	 GSegNode* prev=NULL;
	 while (node && node->start()< r.start) {
		 prev=node;
		 node=node->next;
	 }
	 //prev is the last segment starting before r
	 if (prev) {
		 if (prev->end()>=r.start)  // r overlaps prev segment
			d=r.start - prev->start();
		 if (d==0)
			clearTo(prev); //clear all nodes including prev
	 }

     if (last_pos!=r.start) {
    	 last_pos=r.start;
    	 last_dist=d;
     }
     mergeRead(r);
	 return d;
 }


};


struct RDistanceData {
  GVec<GSegList> fsegs; //forward strand segs for each sample
  GVec<GSegList> rsegs; //reverse strand segs for each sample
  void init(int num_samples) {
	  fsegs.Resize(num_samples);
	  rsegs.Resize(num_samples);
	  this->reset();
  }
  void reset() {
	  for (int i=0;i<fsegs.Count();i++) {
           fsegs[i].reset();
           rsegs[i].reset();
	  }
  }
};

RDistanceData rspacing;

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
class SPData { // Same Point data
    bool settled; //real SP data, owns its r data and deallocates it on destroy
  public:
	int64_t accYC; //if any merged records had YC tag, their values are accumulated here
	int64_t accYX; //if any merged records had YX tag, their values are accumulated here
	int64_t maxYD; //max bundle extent upstream in all samples (0 when no preceding overlapping reads)
	GBitVec* samples; //which samples were the collapsed ones coming from
	                 //number of bits set will be stored as YX:i:(samples.count()+accYX)
	int dupCount; //duplicity count - how many single-alignments were merged into r
	              // will be stored as tag YC:i:(dupCount+accYC)
	GSamRecord* r;
	char tstrand; //'-','+' or '.'
    SPData(GSamRecord* rec=NULL):settled(false), accYC(0), accYX(0), maxYD(0),samples(NULL),
    		dupCount(0), r(rec), tstrand('.') {
    	if (r!=NULL) tstrand=r->spliceStrand();
    }

    ~SPData() {
    	if (settled && r!=NULL) delete r;
    	if (samples!=NULL) delete samples;
    }

    void detach(bool dontFree=true) { settled=!dontFree; }
    //detach(true) must never be called before settle()

    void settle(TInputRecord& trec) { //becomes a standalone SPData record
    	// duplicates the current record
    	settled=true;
    	trec.disown();
    	if (samples==NULL)
    		samples=new GBitVec(inRecords.freaders.Count());

    	if (trec.tbMerged) {
    		accYC=r->tag_int("YC", 1);
    		accYX=r->tag_int("YX", 1);
    		maxYD=r->tag_int("YD", 0);
    	} else {
    		++dupCount;
    		samples->set(trec.fidx);
    	}
    }

    void dupAdd(TInputRecord& trec) { //merge an external SAM record into this one
    	if (!settled) GError("Error: cannot merge a duplicate into a non-settled SP record!\n");
    	GSamRecord& rec=*trec.brec;
    	//WARNING: rec MUST be a "duplicate" of current record r
    	if (trec.tbMerged) {
    		accYC+=rec.tag_int("YC",1);
    		accYX+=rec.tag_int("YX", 1);
    		int64_t vYD=rec.tag_int("YD",0);
    		if (vYD>maxYD) maxYD=vYD; //keep only maximum YD value
    	} else {
    		//avoid collapsing same read alignment duplicated just for pairing reasons
    		if (!samples->test(trec.fidx) || rec.pairOrder()!=r->pairOrder() ||
    				strcmp(r->name(), rec.name())!=0) {
    		   dupCount++;
    		   samples->set(trec.fidx);
    		}
    	}
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
		if (spf!=newspd) { //matches existing SP entry spf, merge
			spf->dupAdd(irec); //update existing SP entry
			delete newspd;
			return;
		} // not a novel SP data
		//else newspd was added as a separate entry of spdlst
	}
	else { // empty list, just add this
		spdlst.Add(newspd);
	}
	newspd->settle(irec); //keep its own SAM record copy
}

void flushPData(GList<SPData>& spdlst){ //write spdata to outfile
  if (spdlst.Count()==0) return;
  // write SAM records in spdata to outfile
  for (int i=0;i<spdlst.Count();++i) {
	  SPData& spd=*(spdlst.Get(i));
	  int64_t accYC=spd.accYC+spd.dupCount;
	  int64_t accYX=spd.accYX;
	  int dSamples=spd.samples->count(); //this only has direct, non-TieBrush samples
	  accYX+=dSamples;
	  if (accYC>1) spd.r->add_int_tag("YC", accYC);
	  if (accYX>1) spd.r->add_int_tag("YX", accYX);
	  int dmax=spd.maxYD;
	  for(int s=spd.samples->find_first();s>=0;
	    		s=spd.samples->find_next(s)) {
	    	if (spd.tstrand=='+' || spd.tstrand=='.') {
	    	   int r=rspacing.fsegs[s].processRead(*spd.r);
	    	   if (r>dmax) dmax=r;
	    	}
	    	if (spd.tstrand=='-' || spd.tstrand=='.') {
	    	   int r=rspacing.rsegs[s].processRead(*spd.r);
	    	   if (r>dmax) dmax=r;
	    	}
	  } //for each bit index/sample
	  spd.maxYD=dmax;
	  if (spd.maxYD>0) spd.r->add_int_tag("YD", spd.maxYD);
	  else spd.r->remove_tag("YD");
	  outfile->write(spd.r);
	  outCounter++;
  }
  spdlst.Clear();
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	int numSamples=inRecords.start();
	if (outfname.is_empty()) outfname="-";
	GSamFileType oftype=(outfname=="-") ?
			GSamFile_SAM : GSamFile_BAM;
	outfile=new GSamWriter(outfname, inRecords.header(), oftype);
	rspacing.init(numSamples);
	TInputRecord* irec=NULL;
	GSamRecord* brec=NULL;

	GList<SPData> spdata(true, true, true); //list of Same Position data, with all possibly merged records
	bool newChr=false;
	int prev_pos=-1;
	int prev_tid=-1;
	while ((irec=inRecords.next())!=NULL) {
		 brec=irec->brec;
		 if (brec->isUnmapped()) continue;
		 inCounter++;
		 int tid=brec->refId();
		 int pos=brec->start; //1-based

		 if (tid!=prev_tid) {
			 if (prev_tid!=-1) newChr=true;
			 prev_tid=tid;
			 prev_pos=-1;
		 }
		 if (pos!=prev_pos) { //new position
			 flushPData(spdata); //also adds read data to rspacing
			 prev_pos=pos;
		 }
		 if (newChr) {
			 rspacing.reset();
			 newChr=false;
		 }
		 addPData(*irec, spdata);
	}
    flushPData(spdata);
	delete outfile;
	inRecords.stop();
    //if (verbose) {
    double p=100.00 - (double)(outCounter*100.00)/(double)inCounter;
    GMessage("%ld input records written as %ld (%.2f%% reduction)\n", inCounter, outCounter, p);
    //}
}
// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
	GArgs args(argc, argv, "help;debug;verbose;version;cigar;clip;exon;CPEDVho:");
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
		inRecords.addFile(ifn);
	}
}


