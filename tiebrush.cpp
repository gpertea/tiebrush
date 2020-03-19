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

//keep track of all SAM alignments starting at the same coordinate
// that were merged into a single alignment
class SPData {
    bool linked; //not an independed copy, do not free r on destroy
	GBitVec samples; //which samples were the collapsed ones coming from
  public:
	GSamRecord* r;
	int dup; //duplicity count - how many alignments were "the same" with r
    SPData(GSamRecord* rec=NULL, bool standalone=false):linked(~standalone),
    		samples(inRecords.readers.Count()), r(NULL), dup(0) {
    	if (linked) r=rec;
    	else {
    		if (rec!=NULL)
    		  r=new GSamRecord(*rec); //FIXME: full copy constructor needed!
    		else { //should never happen
    			r=new GSamRecord();
    			GMessage("Warning: standalone blank SPData object created!\n");
    		}
    	}
    }
    //TODO: add < and == operators accordingly (for merging)
};

void processOptions(int argc, char* argv[]);

GSamWriter* outfile=NULL;
GStr outfname;

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
	 GPVec<SPData> spdata; //list of Same Position data, with all possibly merged records
	 bool more_alns=true;
	 int prev_pos=0;
	 int prev_tid=-1;
	 while ((irec=inRecords.next())!=NULL) {
		 brec=irec->brec;
		 if (brec->isUnmapped()) continue;
		 int tid=brec->refId();
		 int pos=brec->start; //1-based



	 }

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
			GError("Error: only one merge strategy is accepted.\n");
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


