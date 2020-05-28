#include "GSam.h"
#include "GArgs.h"
#include "tmerge.h"
#include "GBitVec.h"

#include <set>
#include <functional>
#include <iostream>
#include <map>
#include <fstream>

#define VERSION "0.0.1"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-o out.bed] in.bam\n"
                  " Other options: \n"
                  "  --min_nh,-N   : minimum NH score (if available) to include when reporting coverage\n"
                  "  --min_qual,-Q    : minimum mapping quality to include when reporting coverage\n";

TInputFiles inRecords;

GStr outfname;

uint64_t inCounter=0;
uint64_t outCounter=0;

bool debugMode=false;
bool verbose=false;

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
            cov_ss<<al_hdr->target_name[tid]<<"\t"<<start<<"\t"<<end<<"\t"<<cov<<std::endl;
        }
        else{
            fprintf(stdout,"%s\t%d\t%d\t%d\n",al_hdr->target_name[tid],start,end,cov);
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

int getYC(bam1_t *in_rec){
    uint8_t* ptr_yc_1=bam_aux_get(in_rec,"YC");
    if(ptr_yc_1){
        return bam_aux2i(ptr_yc_1);
    }
    else{
        return 1;
    }
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    inRecords.setup(VERSION, argc, argv);
    processOptions(argc, argv);
    inRecords.start();
    if (outfname.is_empty()) outfname="-";
    GSamFileType oftype=(outfname=="-") ?
                        GSamFile_SAM : GSamFile_BAM;
    covInterval.set_outfname(outfname.chars());

    TInputRecord* irec=NULL;
    GSamRecord* brec=NULL;
    //bool newChr=false;
    int prev_tid=-1;
    while ((irec=inRecords.next())!=NULL) {
        brec=irec->brec;
        if (brec->isUnmapped()) continue;
        int tid=brec->refId();
        if (tid!=prev_tid) {
            cleanPriorityQueue(inRecords.header());
            prev_tid=tid;
        }
        int accYC = getYC(brec->get_b());
        incCov(inRecords.header(),brec->get_b(),accYC);
    }
    cleanPriorityQueue(inRecords.header());
    covInterval.write(inRecords.header());

    inRecords.stop();
}
// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVho:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
        GMessage(USAGE);
        exit(1);
    }

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
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    outfname=args.getOpt('o');
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        inRecords.Add(ifn);
    }
}