#include "GSam.h"
#include "GArgs.h"
#include "GStr.h"

const char* USAGE="Usage:\n  tiebrush [-o <outfile>.bam] in1.bam in2.bam ...\n";

enum OutType {
  outFASTQ,
  outFASTA,
  outGFF
};

OutType out_type=outFASTQ;
bool all_reads=false;
bool mapped_only=false;

void showfastq(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped() && !all_reads) return;
  if (mapped_only && rec.isUnmapped()) return;
  char* qseq=rec.sequence();
  fprintf(fout, "@%s\n%s\n", rec.name(), qseq);
  GFREE(qseq);
  qseq=rec.qualities();
  fprintf(fout, "+\n%s\n",qseq);
  GFREE(qseq);
}

void showfasta(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped() && !all_reads) return;
  if (mapped_only && rec.isUnmapped()) return;
  char* qseq=rec.sequence();
  char* alt_name=rec.tag_str("XA");
  if (alt_name)
	fprintf(fout, ">%s alt_name=%s\n%s\n", rec.name(), alt_name, qseq);
  else
    fprintf(fout, ">%s\n%s\n", rec.name(), qseq);
  GFREE(qseq);
}

void showgff(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped()) return;
  char tstrand=rec.spliceStrand();
  fprintf(fout, "%s\tbam\tmRNA\t%d\t%d\t.\t%c\t.\tID=%s\n", rec.refName(),
         rec.start, rec.end, tstrand, rec.name());
  for (int i=0;i<rec.exons.Count();i++) {
    fprintf(fout, "%s\tbam\texon\t%d\t%d\t.\t%c\t.\tParent=%s\n", rec.refName(),
           rec.exons[i].start, rec.exons[i].end, tstrand, rec.name());
     }
}

int main(int argc, char *argv[])  {
    GArgs args(argc, argv, "fasta;fastq;gff;all;help;mapped-only;ref="
        "hAMGaqo:r:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
      GMessage(USAGE);
      return 1;
      }
    args.printCmdLine(stderr);

    all_reads=(args.getOpt('A') || args.getOpt("all"));
    mapped_only=(args.getOpt('M') || args.getOpt("mapped-only"));
    if ((all_reads && mapped_only)) {
        GError("Error: incompatible options !\n");
        }

    if (args.getOpt('a') || args.getOpt("fasta"))
       out_type=outFASTA;
    else if (args.getOpt('G') || args.getOpt("gff")) {
       out_type=outGFF;
       mapped_only=true;
       }
    char* cram_ref=NULL;
    cram_ref=args.getOpt('r');
    if (cram_ref==NULL) cram_ref=args.getOpt("ref");
    char* fname=args.nextNonOpt();
    if (fname==NULL || fname[0]==0) {
        GMessage(USAGE);
        return 1;
    }
    GSamReader samreader(fname, cram_ref);
    //GSamReader samreader(fname);
    FILE* fout=stdout;
    char* outfname=args.getOpt('o');
    if (outfname) {
       fout=fopen(outfname, "w");
       if (fout==NULL) {
           fprintf(stderr, "Error creating output file %s\n", outfname);
           return 2;
       }
    }
    GSamRecord *aln=NULL;
    if (out_type==outFASTA) {
        while ((aln=samreader.next())!=NULL) {
           showfasta(*aln, fout);
           delete aln;
        }
    }
    else if (out_type==outGFF) {
        while ((aln=samreader.next())!=NULL) {
           showgff(*aln, fout);
           delete aln;
        }
    }
      else {
        while ((aln=samreader.next())!=NULL) {
          showfastq(*aln, fout);
          delete aln;
        }
    }
    if (fout!=stdout) fclose(fout);
    return 0;
}
