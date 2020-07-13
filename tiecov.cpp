#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "GSam.h"

#define VERSION "0.0.4"

const char* USAGE="TieCov v" VERSION " usage:\n"
" tiecov [-b out.flt.bam] [-c out.bedgraph] [-j out.junctions.bed] in.bam\n"
" Other options: \n"
"  -N   : maximum NH score (if available) to include when reporting coverage\n"
"  -Q   : minimum mapping quality to include when reporting coverage\n";

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
} filters;

GStr covfname, jfname, bfname, infname;
FILE* boutf=NULL;
FILE* coutf=NULL;
FILE* joutf=NULL;

bool debugMode=false;
bool verbose=false;
int juncCount=0;

struct CJunc {
	int start, end;
	char strand;
	uint64_t dupcount;
	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		return (strand==a.strand && start==a.start && end==a.end);
	}
	bool operator<(const CJunc& a) {
		if (start==a.start) return (end<a.end);
		else return (start<a.start);
	}

	void add(CJunc& j) {
       dupcount+=j.dupcount;
	}

	void write(FILE* f, const char* chr) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				chr, start-1, end, juncCount, (long)dupcount, strand);
	}
};

GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount) {
	char strand = r.spliceStrand();
	if (strand!='+' && strand!='-') return;
	for (int i=1;i<r.exons.Count();i++) {
		CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
				dupcount);
		int ei;
		int r=junctions.AddIfNew(j, &ei);
		if (r==-1) {//existing junction, update
			junctions[ei].add(j);
		}
	}
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}

void processOptions(int argc, char* argv[]);


//b_start MUST be passed 1-based
void addCov(GSamRecord& r, int dupCount, GVec<uint64_t>& bcov, int b_start) {
	bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
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

//b_start MUST be passed 1-based
void flushCoverage(sam_hdr_t* hdr, GVec<uint64_t>& bcov,  int tid, int b_start) {
  if (tid<0 || b_start<=0) return;
  int i=0;
  b_start--; //to make it 0-based;
  while (i<bcov.Count()) {
     uint64_t icov=bcov[i];
     int j=i+1;
     while (j<bcov.Count() && icov==bcov[j]) {
    	 j++;
     }
     if (icov!=0)
       fprintf(coutf, "%s\t%d\t%d\t%ld\n", hdr->target_name[tid], b_start+i, b_start+j, (long)icov);
     i=j;
  }
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
	processOptions(argc, argv);
    //htsFile* hts_file=hts_open(infname.chars(), "r");
    //if (hts_file==NULL)
    //   GError("Error: could not open alignment file %s \n",infname.chars());
	GSamReader samreader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
    if (!covfname.is_empty()) {
       if (covfname=="-" || covfname=="stdout")
    	   coutf=stdout;
       else {
          coutf=fopen(covfname.chars(), "w");
          if (coutf==NULL) GError("Error creating file %s\n",
        		  covfname.chars());
          fprintf(coutf, "track type=bedGraph\n");
       }
    }
    if (!jfname.is_empty()) {
       joutf=fopen(jfname.chars(), "w");
       if (joutf==NULL) GError("Error creating file %s\n",
        		  jfname.chars());
       fprintf(joutf, "track name=junctions\n");
    }
    int prev_tid=-1;
    GVec<uint64_t> bcov(2048*1024);
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based
    GSamRecord brec;
	while (samreader.next(brec)) {
		    int nh = brec.tag_int("NH");
		    if(nh>filters.max_nh)  continue;
		    if (brec.mapq()<filters.min_qual) continue;
		    int endpos=brec.end;
		    if (brec.refId()!=prev_tid || (int)brec.start>b_end) {
		    	if (coutf)
			       flushCoverage(samreader.header() , bcov, prev_tid, b_start);
			    if (joutf)
			       flushJuncs(joutf, samreader.refName(prev_tid));
			    b_start=brec.start;
			    b_end=endpos;
			    if (coutf) {
			      bcov.setCount(0);
			      bcov.setCount(b_end-b_start+1);
			    }
			    prev_tid=brec.refId();
		    } else { //extending current bundle
			    if (b_end<endpos) {
				    b_end=endpos;
				    bcov.setCount(b_end-b_start+1, (int)0);
			    }
		    }
		    int accYC = brec.tag_int("YC", 1);
		    if (coutf)
		       addCov(brec, accYC, bcov, b_start);
		    if (joutf && brec.exons.Count()>1) {
		    	addJunction(brec, accYC);
		    }
	} //while GSamRecord emitted
	if (coutf) {
       flushCoverage(samreader.header(), bcov, prev_tid, b_start);
       if (coutf!=stdout) fclose(coutf);
	}
	if (joutf) {
		flushJuncs(joutf, samreader.refName(prev_tid));
		fclose(joutf);
	}

}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVhc:j:b:N:Q:");
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
    covfname=args.getOpt('c');
    jfname=args.getOpt('j');
    bfname=args.getOpt('b');
    if (args.startNonOpt()!=1) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();
}
