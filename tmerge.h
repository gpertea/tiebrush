#ifndef TIEBRUSH_TMERGE_H_
#define TIEBRUSH_TMERGE_H_
#include "GStr.h"
#include "GVec.hh"
#include "GList.hh"
//#include "rlink.h"
#include "GSam.h"
#include "htslib/khash.h"

struct TSamReader {
	GStr fname;
	GSamReader* samreader;
	bool tbMerged; //based on the header, is the file a product of TieBrush?
	TSamReader(const char* fn=NULL, GSamReader* samr=NULL):
		fname(fn), samreader(samr), tbMerged(false) {}
	~TSamReader() {
		delete samreader;
	}
};

struct TInputRecord {
	GSamRecord* brec;
	int fidx; //file index in files and readers
	bool tbMerged; //is it from a TieBrush generated file?
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
    void disown() {
    	brec=NULL;
    }
	TInputRecord(GSamRecord* b=NULL, int i=0, bool tb_merged=false):brec(b),
			fidx(i),tbMerged(tb_merged) {}
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
	GPVec<TSamReader> freaders;
	void addFile(const char* fn);
	bool addSam(GSamReader* r, int fidx); //update mHdr data
	GList<TInputRecord> recs; //next record for each
	TInputFiles():crec(NULL), mHdr(NULL), pg_ver(NULL), pg_args(),
			freaders(true), recs(true, true, true) { }

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

	~TInputFiles() {
		GFREE(pg_ver);
		sam_hdr_destroy(mHdr);
	}

	int count() { return freaders.Count(); }
	int start(); //open all files, load 1 record from each
	TInputRecord* next();
	void stop(); //
};


#endif /* TIEBRUSH_TMERGE_H_ */
