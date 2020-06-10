#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "htslib/sam.h"

#define VERSION "0.0.3"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-o out.bedgraph] in.bam\n"
                  " Other options: \n"
                  "  -N   : maximum NH score (if available) to include when reporting coverage\n"
                  "  -Q   : minimum mapping quality to include when reporting coverage\n";

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
} filters;

GStr outfname, infname;
FILE* outf=NULL;
bool debugMode=false;
bool verbose=false;

void processOptions(int argc, char* argv[]);

void addCov(bam1_t *in_rec, int dupCount, GVec<uint16_t>& bcov, int b_start) {
    int pos=in_rec->core.pos; // 0-based
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
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bcov[pos-b_start]+=dupCount;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}


int getYC(bam1_t *in_rec){
    uint8_t* ptr_yc_1=bam_aux_get(in_rec,"YC");
    if(ptr_yc_1){
        return bam_aux2i(ptr_yc_1);
    }
    else{
        return 1;
    }
}

int getNH(bam1_t *in_rec){
    uint8_t* ptr_nh_1=bam_aux_get(in_rec,"NH");
    if(ptr_nh_1){
        return bam_aux2i(ptr_nh_1);
    }
    else{
        return -1;
    }
}

void flushCoverage(sam_hdr_t* hdr, GVec<uint16_t>& bcov,  int tid, int b_start) {
  if (tid<0 || b_start<0) return;
  int i=0;
  while (i<bcov.Count()) {
     uint16_t icov=bcov[i];
     int j=i+1;
     while (j<bcov.Count() && icov==bcov[j]) {
    	 j++;
     }
     if (icov!=0)
       fprintf(outf, "%s\t%d\t%d\t%d\n", hdr->target_name[tid], b_start+i, b_start+j, icov);
     i=j;
  }
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	processOptions(argc, argv);
    htsFile* hts_file=hts_open(infname.chars(), "r");
    if (hts_file==NULL)
       GError("Error: could not open alignment file %s \n",infname.chars());
    sam_hdr_t* hdr = sam_hdr_read(hts_file);
    int prev_tid=-1;
    bam1_t* b = bam_init1();
    if (outfname.is_empty() || outfname=="-") outf=stdout;
    else {
    	    outf=fopen(outfname.chars(), "w");
    	    if (outf==NULL) GError("Error creating file %s\n", outfname.chars());
    }
    GVec<uint16_t> bcov(4096*1024);
    int b_end=-1; //bundle end, start
    int b_start=-1;
	    while (sam_read1(hts_file, hdr, b) >= 0) {
		    int nh = getNH(b);
		    if(nh>filters.max_nh) { continue; }
		    if (b->core.qual<filters.min_qual) { continue; }
		    int endpos=bam_endpos(b);
		    if (b->core.tid!=prev_tid || b->core.pos>b_end) {
			    flushCoverage(hdr, bcov, prev_tid, b_start);
			    b_start=b->core.pos;
			    b_end=endpos;
			    bcov.setCount(0);
			    bcov.setCount(b_end-b_start+1);
			    prev_tid=b->core.tid;
		    } else { //extending current bundle
			    if (b_end<endpos) {
				    b_end=endpos;
				    bcov.setCount(b_end-b_start+1, (int)0);
			    }
		    }
		    int accYC = getYC(b);
		    addCov(b, accYC, bcov, b_start);
	    }
    flushCoverage(hdr, bcov, prev_tid, b_start);
    if (outf!=stdout) fclose(outf);
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVho:N:Q:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    GStr max_nh_str=args.getOpt('N');
    if (!max_nh_str.is_empty()) {
        filters.max_nh=max_nh_str.asInt();
    }
    GStr min_qual_str=args.getOpt('Q');
    if (!min_qual_str.is_empty()) {
        filters.min_qual=min_qual_str.asInt();
    }

    debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }
    //verbose=(args.getOpt('v')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    outfname=args.getOpt('o');
    if (args.startNonOpt()==0) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();
}
