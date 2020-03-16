#include "tmerge.h"

void TInputFiles::Add(const char* fn) {
	   GStr sfn(fn);
		if (sfn!="-" && fileExists(fn)<2) {
			    GError("Error: input file %s cannot be found!\n",
			            fn);
		}

		files.Add(sfn);
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
		GSamReader* bamreader=new GSamReader(files[i].chars());
		readers.Add(bamreader);
		GSamRecord* brec=bamreader->next();
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

