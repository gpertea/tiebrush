#include "tmerge.h"

void TInputFiles::Add(const char* fn) {
	   GStr sfn(fn);
		if (sfn!="-" && fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}

		files.Add(sfn);
	}

void TInputFiles::addSam(GSamReader& r) {
	kstring_t hd_line = { 0, 0, NULL };
	int res = sam_hdr_find_hd(r.hdr, &hd_line);
	if (res < 0) GError("Error: failed to get @HD line from header");
	//check for SO:coordinate
	kstring_t str = KS_INITIALIZE;
    // Accept unknown, unsorted, or queryname sort order, but error on coordinate sorted.
    if (!sam_hdr_find_tag_hd(r.hdr, "SO", &str) || !str.s || !strcmp(str.s, "coordinate"))
	        GError("Error: %s file not coordinate-sorted!\n");
    ks_free(&str);
	if (!mHdr.haveHD) {
	   if (sam_hdr_add_lines(mHdr.hdr, hd_line.s, hd_line.l) < 0) {
	        GError("Error: failed to add @HD line to new header");
	   }
	   //also add other useful header lines from the first file
	   if (sam_hdr_find_line_id(r.hdr, "PG", NULL, NULL, &str)) {
		   sam_hdr_add_lines(mHdr.hdr, str.s, str.l);
	   }
	   ks_free(&str);
	   sam_hdr_add_pg(mHdr.hdr, "tiebrush",
	                                    "VN", pg_ver, "CL", pg_args);
	   mHdr.haveHD=true;
	}
    //add every single SQ to mHdr.hdr and mHdr.refs
	kstring_t sq_line = { 0, 0, NULL }, sq_sn = { 0, 0, NULL };

	// Fill in the tid part of the translation table, adding new targets
	// to the merged header
	for (int i = 0; i < sam_hdr_nref(r.hdr); ++i) {
		sq_sn.l = 0;
		res = sam_hdr_find_tag_pos(r.hdr, "SQ", i, "SN", &sq_sn);
		if (res < 0)
			GError("Error: failed to get @SQ SN #%d from header\n", i + 1);
		int trans_tid = sam_hdr_name2tid(mHdr.hdr, sq_sn.s);
		if (trans_tid < -1)
			GError("Error: failed to lookup ref\n");

		if (trans_tid < 0) {
			// Append new entry to out_hdr
			sq_line.l = 0;
			res = sam_hdr_find_line_id(r.hdr, "SQ", "SN", sq_sn.s, &sq_line);
			if (res < 0)
				GError("Error: failed to get @SQ SN:%s from header\n", sq_sn.s);

			trans_tid = sam_hdr_nref(mHdr.hdr);

			res = sam_hdr_add_lines(mHdr.hdr, sq_line.s, sq_line.l);
			if (res < 0)
				GError("Error: failed to add @SQ SN:%s to new header\n", sq_sn.s);
		}
		//TODO add to mHdr.refs
		/* --FIXME
		tbl->tid_trans[i] = trans_tid;

		if (tbl->tid_trans[i] > min_tid) {
			min_tid = tbl->tid_trans[i];
		} else {
			tbl->lost_coord_sort = true;
		}
		*/
	}

	ks_free(&sq_line);
	ks_free(&sq_sn);

}

int TInputFiles::start() {
	if (this->files.Count()==1) {
		//special case, if it's only one file it must be a list (usually)
		GStr fname(this->files.First());
		FILE* flst=fopen(this->files.First().chars(),"r");
		if (flst==NULL) GError("Error: could not open input file %s!\n",
				fname.chars());
		GVec<GStr> infiles;
		char* line=NULL;
		int lcap=5000;
		GMALLOC(line, lcap);
		bool firstline=true;
		//bool isalist=true;
		while (fgetline(line,lcap,flst)) {
			GStr s(line);
			s.trim();
			if (s.length()<2 || s[0]=='#') continue; //skip comments/header in the list file, if any
			if (firstline) {
				if (!fileExists(s.chars())) {
					//it must be a GFF
					//this shouldn't normally happen in mergeMode, with one input file
					if (s.count('\t')<8)
						GError("Error: cannot find file '%s' and %s does not look like GFF!\n", s.chars(), fname.chars());
					break;
				}
				firstline=false;
				files.Clear();
			}
			if (!fileExists(s.chars()))
				GError("Error opening transcript file %s !\n",s.chars());
			fname=s;
			files.Add(fname);
		} //for each line in the list file
		GFREE(line);
		fclose(flst);
	}

	//multi-BAM input
	for (int i=0;i<files.Count();++i) {
		GSamReader* samreader=new GSamReader(files[i].chars());
		addSam(*samreader); //merge SAM headers etc.
		GSamRecord* brec=samreader->next();
		if (brec)
		   recs.Add(new TInputRecord(brec, i));
	}
	return readers.Count();
}

GSamRecord* TInputFiles::next() {
	//must free old current record first
	delete crec;
	crec=NULL;
    if (recs.Count()>0) {
    	crec=recs.Pop();//lowest coordinate
    	GSamRecord* rnext=readers[crec->fidx]->next();
    	if (rnext)
    		recs.Add(new TInputRecord(rnext,crec->fidx));
    	return crec->brec;
    }
    else return NULL;
}

void TInputFiles::stop() {
 for (int i=0;i<readers.Count();++i) {
	 readers[i]->bclose();
 }
}

