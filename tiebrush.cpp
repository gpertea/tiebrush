#include "GSam.h"
#include "GArgs.h"
#include "tmerge.h"
#include "GBitVec.h"

#define VERSION "0.0.1"

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

GSamWriter* outfile=NULL;
GStr outfname;

int cmpFull(GSamRecord& a, GSamRecord& b) {
	//-- CIGAR && MD strings
	if (a.b->core.n_cigar!=b.b->core.n_cigar) return ((int)a.b->core.n_cigar - (int)b.b->core.n_cigar);
	int cigar_cmp=0;
	if (a.b->core.n_cigar>0)
		cigar_cmp=memcmp(bam_get_cigar(a.b) , bam_get_cigar(b.b), a.b->core.n_cigar*sizeof(uint32_t) );
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
	if (a.b->core.n_cigar!=b.b->core.n_cigar) return ((int)a.b->core.n_cigar - (int)b.b->core.n_cigar);
	if (a.b->core.n_cigar==0) return 0;
	return memcmp(bam_get_cigar(a.b) , bam_get_cigar(b.b), a.b->core.n_cigar*sizeof(uint32_t) );
}

int cmpCigarClip(GSamRecord& a, GSamRecord& b) {
	uint32_t a_clen=a.b->core.n_cigar;
	uint32_t b_clen=b.b->core.n_cigar;
	uint32_t* a_cstart=bam_get_cigar(a.b);
	uint32_t* b_cstart=bam_get_cigar(b.b);
	while (a_clen>0 &&
			((*a_cstart) && BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { a_cstart++; a_clen--; }
	while (a_clen>0 &&
			(a_cstart[a_clen-1] && BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) a_clen--;
	while (b_clen>0 &&
			((*b_cstart) && BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { b_cstart++; b_clen--; }
	while (b_clen>0 &&
			(b_cstart[b_clen-1] && BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) b_clen--;
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
    bool detached; //if detached, r is not deleted on destroy
	GBitVec samples; //which samples were the collapsed ones coming from
  public:
	GSamRecord* r;
	int dupCount; //duplicity count - how many alignments were "the same" with r
    SPData(GSamRecord* rec=NULL, bool dupRec=false):detached(false),
    		samples(inRecords.readers.Count()), r(NULL), dupCount(0) {
    	if (dupRec) {
    		if (rec!=NULL)
    		  r=new GSamRecord(*rec); //copy constructor
    		else  //should never happen
    		  GError("Error: record duplication of NULL record\n");
    	}
    	else r=rec;
    }
    ~SPData() {
    	if (!detached && r) delete r;
    }
    void detach(bool dontFree=true) { detached=dontFree; }

    bool operator<(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId()) return (r->refId()<b.r->refId());
    	//NOTE: already assuming that start&end must match, no matter the merge strategy
    	if (r->start!=b.r->start) return (r->start<b.r->start);
    	if (r->end!=b.r->end) return (r->end<b.r->end);

    	if (mrgStrategy==tMrgStratFull) return (cmpFull(*r, *(b->r))<0);
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b->r))<0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b->r))<0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b->r))<0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }

    bool operator==(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId() || r->start!=b.r->start ||
    			r->end!=b.r->end) return false;
    	if (mrgStrategy==tMrgStratFull) return (cmpFull(*r, *(b->r))==0);
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b->r))==0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b->r))==0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b->r))==0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }
};

void processOptions(int argc, char* argv[]);

void addPData(TInputRecord& irec, GList<SPData>& spdata) {
  //add and collapse
}

void flushPData(GList<SPData>& spdata){ //write spdata to outfile
  //TODO: write SAM records in spdata to outfile
  spdata.Clear();
}


bool debugMode=false;
bool verbose=false;

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	inRecords.start();
	GSamFileType oftype=(outfname.is_empty() || outfname=="-") ?
			GSamFile_SAM : GSamFile_BAM;
	outfile=new GSamWriter(outfname, inRecords.header(), oftype);

	 TInputRecord* irec=NULL;
	 GSamRecord* brec=NULL;
	 GList<SPData> spdata(true, true, true); //list of Same Position data, with all possibly merged records
	 bool newChr=false;  //flush spdata
	 bool newPos=false; //flush spdata
	 int prev_pos=-1;
	 int prev_tid=-1;
	 while ((irec=inRecords.next())!=NULL) {
		 brec=irec->brec;
		 if (brec->isUnmapped()) continue;
		 int tid=brec->refId();
		 int pos=brec->start; //1-based
		 if (tid!=prev_tid) {
			 prev_tid=tid;
			 newChr=true;
			 prev_pos=-1;
		 }
		 if (pos!=prev_pos) {
			 flushPData(spdata);
			 prev_pos=pos;
		 }
		 addPData(*irec, spdata);
	 }
     flushPData(spdata);
	delete outfile;
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
	GStr outfname=args.getOpt('o');
}


