#ifndef TIEBRUSH_TMERGE_H_
#define TIEBRUSH_TMERGE_H_
#include "GStr.h"
#include "GVec.hh"
#include "GList.hh"
//#include "rlink.h"
#include "GSam.h"
#include "htslib/khash.h"

struct TInputRecord {
	GSamRecord* brec;
	int fidx; //index in files and readers
	bool operator<(TInputRecord& o) {
		 //decreasing location sort
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 int r1_tid=r1.refId();
		 int r2_tid=r2.refId();
		 //int refcmp=strcmp(r1.refName(),r2.refName()); //should use the refID comparison instead
		 if (r1_tid==r2_tid) {
		 //higher coords first
			if (r1.start!=r2.start)
				 return (r1.start>r2.start);
			else {
				if (r1.end!=r2.end)
				   return (r1.end>r2.end);
				else if (fidx==o.fidx)
						return strcmp(r1.name(), r2.name())>0;
					else return fidx>o.fidx;
			}
		 }
		 else {
			 return (r1_tid>r2_tid);
		 }
	}
	bool operator==(TInputRecord& o) {
		 GSamRecord& r1=*brec;
		 GSamRecord& r2=*(o.brec);
		 return ( r1.refId()==r2.refId() && r1.start==r2.start && r1.end==r2.end
				 && fidx==o.fidx && r1.get_b()->l_data==r2.get_b()->l_data &&
				 memcmp(r1.get_b()->data, r1.get_b()->data, r1.get_b()->l_data)==0
				 );
	}

	TInputRecord(GSamRecord* b=NULL, int i=0):brec(b),fidx(i) {}
	~TInputRecord() {
		delete brec;
	}
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	// use that to check if each input SAM file has the refseqs sorted by coordinate
	sam_hdr_t* mHdr; //merged output header data
	char* pg_ver;
	GStr pg_args;
 public:
	GVec<GStr> files; //same order as readers
	GPVec<GSamReader> readers;
	void addSam(GSamReader* r); //update mHdr data
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), mHdr(NULL), pg_ver(NULL), pg_args(),
			files(), readers(true), recs(true, true, true) { }

	sam_hdr_t* header() { return mHdr; }

	void setup(const char* ver, int argc, char** argv) {
		if (ver) pg_ver=Gstrdup(ver);
		if (argc>0 && argv!=NULL) {
			 for (int i=0;i<argc;i++) {
			   pg_args.append(argv[i]);
			   if (i<argc-1) pg_args.append(' ');
			 }
		}
	}

	~TInputFiles() { GFREE(pg_ver); sam_hdr_destroy(mHdr); }
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	TInputRecord* next();
	void stop(); //
};


#endif /* TIEBRUSH_TMERGE_H_ */
