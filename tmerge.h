#ifndef TIEBRUSH_TMERGE_H_
#define TIEBRUSH_TMERGE_H_
#include "GStr.h"
#include "GList.hh"
//#include "rlink.h"
#include "GSam.h"
#include "htslib/khash.h";

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

struct TMrgSQData {
	char* name;
	int sqlen;
	int tid; //newly assigned tid in the merged output header
    TMrgSQData(const char* n=NULL, int sql=0, int id=0):
    	name(NULL),sqlen(sql),tid(id) {
    	if (n) name=Gstrdup(n);
    }
    ~TMrgSQData() { GFREE(name); }
};

struct TMrgHeader {
    sam_hdr_t    *hdr; //merged header for output
    GVec<TMrgSQData> refs; //collected ref seq data
    bool haveHD;
    //khash_t(c2i) *sq_tids;
    TMrgHeader():hdr(NULL), refs(), haveHD(false) {
    	hdr=sam_hdr_init();
    	if (!hdr) GError("Error: could not create SAM header!\n");
    }
    ~TMrgHeader() {
    	if (!hdr) sam_hdr_destroy(hdr);
    }
};

struct TInputFiles {
 protected:
	TInputRecord* crec;
	// use that to check if each input SAM file has the refseqs sorted by coordinate
	TMrgHeader mHdr; //merged output header data
	char* pg_ver;
	char* pg_args;
 public:
	GVec<GStr> files; //same order as readers
	GPVec<GSamReader> readers;
	void addSam(GSamReader& r); //update mHdr data
	GList<TInputRecord> recs; //next record for each
	TInputFiles(const char* ver, const char* args):crec(NULL), mHdr(), pg_ver(NULL), pg_args(NULL),
			files(), readers(true), recs(true, true, true) {
		if (ver) pg_ver=Gstrdup(ver);
		if (args) pg_args=Gstrdup(args);
	}
	~TInputFiles() { GFREE(pg_ver); GFREE(pg_args); }
	void Add(const char* fn);
	int count() { return files.Count(); }
	int start(); //open all files, load 1 record from each
	GSamRecord* next();
	void stop(); //
};


#endif /* TIEBRUSH_TMERGE_H_ */
