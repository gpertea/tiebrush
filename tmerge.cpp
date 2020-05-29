#include "tmerge.h"

void TInputFiles::Add(const char* fn) {
	   GStr sfn(fn);
		if (sfn!="-" && fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}

		files.Add(sfn);
	}

void TInputFiles::addSam(GSamReader* r) {
	//requirement: all files must have the same number of SQ entries in the same order!
	kstring_t hd_line = { 0, 0, NULL };
	int res = sam_hdr_find_hd(r->header(), &hd_line);
	if (res < 0) GError("Error: failed to get @HD line from header");
	//check for SO:coordinate
	kstring_t str = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(r->header(), "SO", &str)
    		|| !str.s
			|| strcmp(str.s, "coordinate") )
	        GError("Error: %s file not coordinate-sorted!\n", r->fileName());
    ks_free(&hd_line);
	if (mHdr==NULL) { //first file
		mHdr=sam_hdr_dup(r->header());
		sam_hdr_add_pg(mHdr, "TieBrush",
				"VN", pg_ver, "CL", pg_args.chars(), NULL);
		// sam_hdr_rebuild(mHdr); -- is this really needed?
	}
	else { //check if this file has the same SQ entries in the same order
		//if it has more seqs, make it the main header
		int r_numrefs=sam_hdr_nref(r->header());
		int m_nrefs=sam_hdr_nref(mHdr);
		bool swapHdr=(r_numrefs>m_nrefs);
		sam_hdr_t *loHdr = swapHdr ? mHdr : r->header();
		sam_hdr_t *hiHdr = swapHdr ? r->header() : mHdr;
		int loNum = swapHdr ? m_nrefs : r_numrefs;
		//if (r_numrefs!=m_nrefs)
		// GError("Error: file %s has different number of reference sequences (%d)!\n", r->fileName(), r_numrefs);
		for (int i = 0; i < loNum; ++i) {
			str.l = 0;
			res = sam_hdr_find_tag_pos(loHdr, "SQ", i, "SN", &str);
			if (res < 0)
				GError("Error: failed to get @SQ SN #%d from header\n", i + 1);
			int m_tid = sam_hdr_name2tid(hiHdr, str.s);
			if (m_tid < -1)
				GError("Error: unexpected ref lookup failure (%s)!\n", str.s);
			if (m_tid < 0)
				GError("Error: ref %s not seen before!\n", str.s);
			int r_tid = sam_hdr_name2tid(loHdr, str.s);
			if (r_tid != m_tid)
					GError("Error: ref %s from file %s does not have the expected id#!", str.s, r->fileName());
		}
		if (swapHdr) {
			sam_hdr_destroy(mHdr);
			mHdr=sam_hdr_dup(r->header());
		}
	}

	ks_free(&str);
    readers.Add(r);
}

int TInputFiles::start() {
	if (this->files.Count()==1) {
		//special case, if it's only one file it might be a list of file paths
		GStr fname(this->files.First());
		//try to open it as a SAM/BAM/CRAM
		htsFile* hf=hts_open(fname.chars(), "r");
		if (hf==NULL || hf->format.category!=sequence_data) {
			//must be a list
			if (hf) hts_close(hf);
#ifdef _DEBUG
			GMessage("Warning: try to open file %s as a list.\n", fname.chars());
#endif
			FILE* flst=fopen(fname.chars(),"r");
			if (flst==NULL) GError("Error: could not open input file %s!\n",
					fname.chars());
			char* line=NULL;
			int lcap=5000;
			GMALLOC(line, lcap);
			files.Clear();
			while (fgetline(line,lcap,flst)) {
				GStr s(line);
				s.trim();
				if (s.length()<2 || s[0]=='#') continue; //skip comments/header in the list file, if any
				if (!fileExists(s.chars()))
					GError("Error: cannot find alignment file %s !\n",s.chars());
				files.Add(s);
			} //for each line in the list file
			GFREE(line);
			fclose(flst);
		}
		else { // single alignment file
			hts_close(hf);
		}
	}

	for (int i=0;i<files.Count();++i) {
		GSamReader* samreader=new GSamReader(files[i].chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
		addSam(samreader); //merge SAM headers etc.
		GSamRecord* brec=samreader->next();
		if (brec)
		   recs.Add(new TInputRecord(brec, i));
	}
	return readers.Count();
}

TInputRecord* TInputFiles::next() {
	//must free old current record first
	delete crec;
	crec=NULL;
    if (recs.Count()>0) {
    	crec=recs.Pop();//lowest coordinate
    	GSamRecord* rnext=readers[crec->fidx]->next();
    	if (rnext)
    		recs.Add(new TInputRecord(rnext,crec->fidx));
    	//return crec->brec;
    	return crec;
    }
    else return NULL;
}

void TInputFiles::stop() {
 for (int i=0;i<readers.Count();++i) {
	 readers[i]->bclose();
 }
}

