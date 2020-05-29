#include "GArgs.h"
#include "GStr.h"
#include "GVec.hh"
#include "htslib/sam.h"
#include <set>
#include <functional>
#include <iostream>
#include <map>
#include <fstream>

#define VERSION "0.0.2"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-o out.bedgraph] in.bam\n"
                  " Other options: \n"
                  "  -C   : use alternate (faster) coverage calculation\n"
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
bool alesCov=true;

void processOptions(int argc, char* argv[]);

std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> > cov_pq; // since RB+ Tree is sorted by std::greater - smallest position (pair of seqid and pos) will always be at the top (begin())
std::pair<std::map<std::pair<uint16_t,uint32_t>,uint16_t,std::less<std::pair<uint16_t,uint32_t>> >::iterator,bool>  cpq_it;

struct CovInterval{
    int tid=-1;
    int start=-1; // start position
    int end=-1; // end position
    int cov=-1; // current coverage for the interval

    CovInterval() = default;
    ~CovInterval(){
        if(cov_ss.is_open()){
            this->cov_ss.close();
        }
    }

    void set_outfname(std::string out_fname){
        if(out_fname=="-" || out_fname.empty()){
            return;
        }
        this->cov_ss.open(out_fname,std::ios::out);
    }

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
        if(cov_ss.is_open()){
            cov_ss<<al_hdr->target_name[tid]<<"\t"<<start-1<<"\t"<<end<<"\t"<<cov<<std::endl;
        }
        else{
            fprintf(stdout,"%s\t%d\t%d\t%d\n",al_hdr->target_name[tid],start-1,end,cov);
        }
    }
private:
    std::fstream cov_ss;
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
                //pos+=oplen; //WRONG
                break;
            case BAM_CMATCH: // base match - increment coverage
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
    if (alesCov) {
    	covInterval.set_outfname(outfname.chars());
		while (sam_read1(hts_file, hdr, b) >= 0) {
			int nh = getNH(b);
			if(nh>filters.max_nh){continue;}
			if (b->core.qual<filters.min_qual) { continue; }
			int tid=b->core.tid;
			if (tid!=prev_tid) {
				cleanPriorityQueue(hdr);
				prev_tid=tid;
			}
			int accYC = getYC(b);
			incCov(hdr, b, accYC);
		}
		cleanPriorityQueue(hdr);
		covInterval.write(hdr);
    }
    else {
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
				bcov.setCount(b_end-b_start+1, (int)0);
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
    }
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;CDVho:N:Q:");
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
    if (args.getOpt('C'))
       alesCov=false ;
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
