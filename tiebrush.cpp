#include "GSam.h"
#include "GArgs.h"
#include "tmerge.h"

#define VERSION "0.0.1"

const char* USAGE="TieBrush v" VERSION " usage:\n"
" tiebrush [-o <outfile>.bam] list.txt | in1.bam in2.bam ...\n";

void processOptions(GArgs& args);

GSamWriter* outfile=NULL;
GStr outfname;

bool debugMode=false;
bool verbose=false;
TInputFiles bamreader;

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	bamreader.setup(VERSION, argc, argv);
	GArgs args(argc, argv, "help;debug;verbose;version;DVho:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
		GMessage(USAGE);
		return 1;
	}
	bamreader.start();
	GSamFileType oftype=(outfname.is_empty() || outfname=="-") ?
			GSamFile_SAM : GSamFile_BAM;
	outfile=new GSamWriter(outfname, bamreader.header(), oftype);


	delete outfile;
}
// <------------------ main() end -----

void processOptions(GArgs& args) {
	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
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


