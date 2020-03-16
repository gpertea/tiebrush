#ifndef TIEBRUSH_TMERGE_H_
#define TIEBRUSH_TMERGE_H_
#include "GStr.h"
#include "GList.hh"
//#include "rlink.h"
#include "GSam.h"

struct TInputRecord {
	GSamRecord* brec;
	int fidx; //index in files and readers
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 int refcmp=strcmp(r1.refName(),r2.refName()); //FIXME no, should use the refID comparison instead
		 if (refcmp==0) {
		 //higher coords first
			if (r1.start!=r2.start)
				 return (r1.start>r2.start);
			else {
				if (r1.end!=r2.end)
				   return (r1.end>r2.end);
				else if (fidx==o.fidx)
						return (strcmp(r1.name(), r2.name())>0); //FIXME
					else return fidx>o.fidx;
			}
		 }
		 else { //use lexicographic order of ref seqs
			 return (refcmp>0);
		 }
	}
	bool operator==(TInputRecord& o) {
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 //TODO: should use refID comparison instead of refName!
		 return ( strcmp(r1.refName(),r2.refName())==0 && r1.start==r2.start && r1.end==r2.end
				 && fidx==o.fidx && strcmp(r1.name(),r2.name())==0);
	}

	TInputRecord(GSamRecord* b=NULL, int i=0):brec(b),fidx(i) {}
	~TInputRecord() {
		delete brec;
	}
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	// TODO: populated a refseq list/ID from the 1st SAM header and
	// use that to check if each input SAM file has the refseqs sorted
	// in the same order!
 public:
	GVec<GStr> files; //same order as readers
	GPVec<GSamReader> readers;
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), files(), readers(true),
			recs(true, true, true) { }
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	GSamRecord* next();
	void stop(); //
};


#endif /* TIEBRUSH_TMERGE_H_ */
