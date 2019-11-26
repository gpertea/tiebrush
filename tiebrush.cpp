#include "alnmerge.h"

#define VERSION "0.0.5"

#define USAGE "TieBrush v" VERSION " usage:\n\
 stringtie <input.bam ...> [-o <output>.sbam]\n\
"

//---- globals

FILE* f_out=NULL;
GStr outfname;

bool debugMode=false;
bool verbose=false;


TInputFiles bamreader;

int main(int argc, char* argv[]) {

 // == Process arguments.
 GArgs args(argc, argv,
   "debug;help;version;o:");
 args.printError(USAGE, true);

 processOptions(args);

 GVec<GRefData> refguides; // plain vector with data for each chromosome

 int bamcount=bamreader.start(); //setup and open input files
 if (bamcount<1) {
	 GError("%sError: no input files provided!\n",USAGE);
 }


const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";

 // --- here we do the input processing
 gseqNames=GffObj::names; //might have been populated already by gff data
 gffnames_ref(gseqNames);  //initialize the names collection if not guided


 GHash<int> hashread;      //read_name:pos:hit_index => readlist index

 int currentstart=0, currentend=0;
 GStr lastref;
 bool no_ref_used=true;
 int lastref_id=-1; //last seen gseq_id
 // int ncluster=0; used it for debug purposes only

 BundleData bundles[1];
 BundleData* bundle = &(bundles[0]);

 GBamRecord* brec=NULL;
 bool more_alns=true;
 int prev_pos=0;
 bool skipGseq=false;
 while (more_alns) {
	 bool chr_changed=false;
	 int pos=0;
	 const char* refseqName=NULL;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 //delete brec;
	 if ((brec=bamreader.next())!=NULL) {
		 if (brec->isUnmapped()) continue;
		 if (brec->start<1 || brec->mapped_len<10) {
			 if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
					 brec->name(), brec->start, brec->mapped_len);
			 continue;
		 }

		 refseqName=brec->refName();
		 xstrand=brec->spliceStrand(); // tagged strand gets priority
		 if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
			 if(brec->isPaired()) { // read is paired
				 if(brec->pairOrder()==1) { // first read in pair
					 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
					 else xstrand='-';
				 }
				 else {
					 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
					 else xstrand='+';
				 }
			 }
			 else {
				 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
				 else xstrand='-';
			 }
		 }

		 /*
		 if (xstrand=='.' && brec->exons.Count()>1) {
			 no_xs++;
			 continue; //skip spliced alignments lacking XS tag (e.g. HISAT alignments)
		 }
		 // I might still infer strand later */

		 if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
		 pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
		 chr_changed=(lastref.is_empty() || lastref!=refseqName);
		 if (chr_changed) {
			 gseq_id=gseqNames->gseqs.addName(refseqName);
			 prev_pos=0;
		 }
		 if (pos<prev_pos) GError(ERR_BAM_SORT);
		 prev_pos=pos;
		 alncounts[gseq_id]++;
		 nh=brec->tag_int("NH");
		 if (nh==0) nh=1;
		 hi=brec->tag_int("HI");
		 //if (mergeMode) {
		    tinfo=new TAlnInfo(brec->name(), brec->tag_int("ZF"));
		    GStr score(brec->tag_str("ZS"));
		    if (!score.is_empty()) {
		      GStr srest=score.split('|');
		      if (!score.is_empty())
		         tinfo->cov=score.asDouble();
		      score=srest.split('|');
		      if (!srest.is_empty())
		    	 tinfo->fpkm=srest.asDouble();
		      srest=score.split('|');
		      if (!score.is_empty())
		         tinfo->tpm=score.asDouble();
		    }
		//}

		 if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
			 new_bundle=true;
		 }
	 }
	 else { //no more alignments
		 more_alns=false;
		 new_bundle=true; //fake a new start (end of last bundle)
	 }

	 if (new_bundle || chr_changed) {
		 //bgeneids.Clear();
		 hashread.Clear();
		 if (bundle->readlist.Count()>0) { // process reads in previous bundle
			// (readthr, junctionthr, mintranscriptlen are globals)
			bundle->getReady(currentstart, currentend);
            /*
			if (gfasta!=NULL) { //genomic sequence data requested
				GFaSeqGet* faseq=gfasta->fetch(bundle->refseq.chars());
				if (faseq==NULL) {
					GError("Error: could not retrieve sequence data for %s!\n", bundle->refseq.chars());
				}
				bundle->gseq=faseq->copyRange(bundle->start, bundle->end, false, true);
			}
            */
			//Num_Fragments+=bundle->num_fragments;
			//Frag_Len+=bundle->frag_len;
			processBundle(bundle);
			// ncluster++; used it for debug purposes only
		 } //have alignments to process
		 else { //no read alignments in this bundle?
			bundle->Clear();
		 } //nothing to do with this bundle

		 if (chr_changed) {
			 lastref=refseqName;
			 lastref_id=gseq_id;
			 currentend=0;
		 }

		 if (!more_alns) {
				if (verbose) {
					if (Num_Fragments) {
					   printTime(stderr);
					   GMessage(" %g aligned fragments found.\n", Num_Fragments);
					}
					//GMessage(" Done reading alignments.\n");
				}
			 noMoreBundles();
			 break;
		 }
		currentstart=pos;
		currentend=brec->end;
		bundle->refseq=lastref;
		bundle->start=currentstart;
		bundle->end=currentend;
	 } //<---- new bundle started

	 if (currentend<(int)brec->end) {
		 //current read extends the bundle
		 //this might not happen if a longer guide had already been added to the bundle
		 currentend=brec->end;
		 if (guides) { //add any newly overlapping guides to bundle
			 bool cend_changed;
			 do {
				 cend_changed=false;
				 while (ng_end+1<ng && (int)(*guides)[ng_end+1]->start<=currentend) {
					 ++ng_end;
					 //more transcripts overlapping this bundle?
					 if ((int)(*guides)[ng_end]->end>=currentstart) {
						 //it should really overlap the bundle
						 bundle->keepGuide((*guides)[ng_end],
								  &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
						 if(currentend<(int)(*guides)[ng_end]->end) {
							 currentend=(*guides)[ng_end]->end;
							 cend_changed=true;
						 }
				 	 }
				 }
			 } while (cend_changed);
		 }
	 } //adjusted currentend and checked for overlapping reference transcripts
	 GReadAlnData alndata(brec, 0, nh, hi, tinfo);
     bool ovlpguide=bundle->evalReadAln(alndata, xstrand);
     }
 } //for each read alignment

 //cleaning up
 delete brec;
 //bamreader.bclose();
 bamreader.stop(); //close all BAM files

 if (verbose) {
    printTime(stderr);
    GMessage(" Done.\n");
 }

#ifdef B_DEBUG
 fclose(dbg_out);
#endif

 fclose(f_out);
 if (c_out && c_out!=stdout) fclose(c_out);

 if (!keepTempFiles) {
   tmp_path.chomp('/');
   remove(tmp_path);
 }


 gffnames_unref(gseqNames); //deallocate names collection


} // -- END main

//----------------------------------------
char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

void processOptions(GArgs& args) {


	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}

	 longreads=(args.getOpt('L')!=NULL);
	 if(longreads) {
		 bundledist=0;
		 singlethr=1.5;
	 }


	if (args.getOpt("conservative")) {
	  isofrac=0.05;
	  singlethr=4.75;
	  readthr=1.5;
	  trim=false;
	}

	if (args.getOpt('t')) {
		trim=false;
	}

	if (args.getOpt("fr")) {
		fr_strand=true;
	}
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}

	 debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	 forceBAM=(args.getOpt("bam")!=NULL); //assume the stdin stream is BAM instead of text SAM
	 mergeMode=(args.getOpt("merge")!=NULL);
	 keepTempFiles=(args.getOpt("keeptmp")!=NULL);
	 //adaptive=!(args.getOpt('d')!=NULL);
	 verbose=(args.getOpt('v')!=NULL);
	 if (verbose) {
	     fprintf(stderr, "Running StringTie " VERSION ". Command line:\n");
	     args.printCmdLine(stderr);
	 }
	 //complete=!(args.getOpt('i')!=NULL);
	 // trim=!(args.getOpt('t')!=NULL);
	 includesource=!(args.getOpt('z')!=NULL);
	 //EM=(args.getOpt('y')!=NULL);
	 //weight=(args.getOpt('w')!=NULL);

	 GStr s=args.getOpt('m');
	 if (!s.is_empty()) {
	   mintranscriptlen=s.asInt();
	   if (!mergeMode) {
		   if (mintranscriptlen<30)
			   GError("Error: invalid -m value, must be >=30)\n");
	   }
	   else if (mintranscriptlen<0) GError("Error: invalid -m value, must be >=0)\n");
	 }
	 else if(mergeMode) mintranscriptlen=50;

	 /*
	 if (args.getOpt('S')) {
		 // sensitivitylevel=2; no longer supported from version 1.0.3
		 sensitivitylevel=1;
	 }
	*/

	 s=args.getOpt("rseq");
	 if (s.is_empty())
		 s=args.getOpt('S');
	 if (!s.is_empty()) {
		 gfasta=new GFastaDb(s.chars());
	 }

     s=args.getOpt('x');
     if (!s.is_empty()) {
    	 //split by comma and populate excludeGSeqs
    	 s.startTokenize(" ,\t");
    	 GStr chrname;
    	 while (s.nextToken(chrname)) {
    		 excludeGseqs.Add(chrname.chars(),new int(0));
    	 }
     }

     /*
	 s=args.getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 GMessage("sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(sensitivitylevel>3) {
			 sensitivitylevel=3;
			 GMessage("sensitivity level out of range: setting sensitivity level at 2\n");
		 }
	 }
	*/


	 s=args.getOpt('p');
	 if (!s.is_empty()) {
		   num_cpus=s.asInt();
		   if (num_cpus<=0) num_cpus=1;
	 }

	 s=args.getOpt('a');
	 if (!s.is_empty()) {
		 junctionsupport=(uint)s.asInt();
	 }

	 s=args.getOpt('j');
	 if (!s.is_empty()) junctionthr=s.asInt();

	 s=args.getOpt('c');
	 if (!s.is_empty()) {
		 readthr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -c value, must be >=0.001)\n");
		 }
	 }
	 else if(mergeMode) readthr=0;



	 s=args.getOpt('g');
	 if (!s.is_empty()) {
		 bundledist=s.asInt();
		 if(bundledist>runoffdist) runoffdist=bundledist;
	 }
	 else if(mergeMode) bundledist=250; // should figure out here a reasonable parameter for merge

	 s=args.getOpt('F');
	 if (!s.is_empty()) {
		 fpkm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) fpkm_thr=0;

	 s=args.getOpt('T');
	 if (!s.is_empty()) {
		 tpm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) tpm_thr=0;

	 s=args.getOpt('l');
	 if (!s.is_empty()) label=s;
	 else if(mergeMode) label="MSTRG";

	 s=args.getOpt('f');
	 if (!s.is_empty()) {
		 isofrac=(float)s.asDouble();
		 if(isofrac>=1) GError("Miminum isoform fraction (-f coefficient: %f) needs to be less than 1\n",isofrac);
	 }
	 else if(mergeMode) isofrac=0.01;
	 s=args.getOpt('M');
	 if (!s.is_empty()) {
		 mcov=(float)s.asDouble();
	 }

	 genefname=args.getOpt('A');
	 if(!genefname.is_empty()) {
		 geneabundance=true;
	 }

	 tmpfname=args.getOpt('o');

	 // coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
	 s=args.getOpt('s');
	 if (!s.is_empty()) {
		 singlethr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -s value, must be >=0.001)\n");
		 }
	 }

	 if (args.getOpt('G')) {
	   guidegff=args.getOpt('G');
	   if (fileExists(guidegff.chars())>1)
	        guided=true;
	   else GError("Error: reference annotation file (%s) not found.\n",
	             guidegff.chars());
	 }

	 enableNames=(args.getOpt('E')!=NULL);

	 retained_intron=(args.getOpt('i')!=NULL);

	 nomulti=(args.getOpt('u')!=NULL);

	 //isunitig=(args.getOpt('U')!=NULL);

	 eonly=(args.getOpt('e')!=NULL);
	 if(eonly && mergeMode) {
		 eonly=false;
		 includecov=true;
	 }
	 else if(eonly && !guided)
		 GError("Error: invalid -e usage, GFF reference not given (-G option required).\n");


	 ballgown_dir=args.getOpt('b');
	 ballgown=(args.getOpt('B')!=NULL);
	 if (ballgown && !ballgown_dir.is_empty()) {
		 GError("Error: please use either -B or -b <path> options, not both.");
	 }
	 if ((ballgown || !ballgown_dir.is_empty()) && !guided)
		 GError("Error: invalid -B/-b usage, GFF reference not given (-G option required).\n");

	 /* s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else { */
		 s=args.getOpt('C');
		 if (!s.is_empty()) {
			 if(!guided) GError("Error: invalid -C usage, GFF reference not given (-G option required).\n");
			 c_out=fopen(s.chars(), "w");
			 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 }
	 //}
	int numbam=args.startNonOpt();
#ifndef GFF_DEBUG
	if (numbam < 1 ) {
	 	 GMessage("%s\nError: no input file provided!\n",USAGE);
	 	 exit(1);
	}
#endif
	const char* ifn=NULL;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		bamreader.Add(ifn);
	}
	//deferred creation of output path
	outfname="stdout";
	out_dir="./";
	 if (!tmpfname.is_empty() && tmpfname!="-") {
		 if (tmpfname[0]=='.' && tmpfname[1]=='/')
			 tmpfname.cut(0,2);
		 outfname=tmpfname;
		 int pidx=outfname.rindex('/');
		 if (pidx>=0) {//path given
			 out_dir=outfname.substr(0,pidx+1);
			 tmpfname=outfname.substr(pidx+1);
		 }
	 }
	 else { // stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname.tr(":","-");
		tmpfname+='.';
		tmpfname+=stime;
	 }
	 if (out_dir!="./") {
		 if (fileExists(out_dir.chars())==0) {
			//directory does not exist, create it
			if (Gmkdir(out_dir.chars()) && !fileExists(out_dir.chars())) {
				GError("Error: cannot create directory %s!\n", out_dir.chars());
			}
		 }
	 }
	 if (!genefname.is_empty()) {
		 if (genefname[0]=='.' && genefname[1]=='/')
		 			 genefname.cut(0,2);
	 //attempt to create the gene abundance path
		 GStr genefdir("./");
		 int pidx=genefname.rindex('/');
		 if (pidx>=0) { //get the path part
				 genefdir=genefname.substr(0,pidx+1);
				 //genefname=genefname.substr(pidx+1);
		 }
		 if (genefdir!="./") {
			 if (fileExists(genefdir.chars())==0) {
				//directory does not exist, create it
				if (Gmkdir(genefdir.chars()) || !fileExists(genefdir.chars())) {
					GError("Error: cannot create directory %s!\n", genefdir.chars());
				}
			 }
		 }

	 }

	 { //prepare temp path
		 GStr stempl(out_dir);
		 stempl.chomp('/');
		 stempl+="/tmp.XXXXXXXX";
		 char* ctempl=Gstrdup(stempl.chars());
	     Gmktempdir(ctempl);
	     tmp_path=ctempl;
	     tmp_path+='/';
	     GFREE(ctempl);
	 }

	 tmpfname=tmp_path+tmpfname;

	 if (ballgown) ballgown_dir=out_dir;
	   else if (!ballgown_dir.is_empty()) {
			ballgown=true;
			ballgown_dir.chomp('/');ballgown_dir+='/';
			if (fileExists(ballgown_dir.chars())==0) {
				//directory does not exist, create it
				if (Gmkdir(ballgown_dir.chars()) && !fileExists(ballgown_dir.chars())) {
					GError("Error: cannot create directory %s!\n", ballgown_dir.chars());
				}

			}
	   }
#ifdef B_DEBUG
	 GStr dbgfname(tmpfname);
	 dbgfname+=".dbg";
	 dbg_out=fopen(dbgfname.chars(), "w");
	 if (dbg_out==NULL) GError("Error creating debug output file %s\n", dbgfname.chars());
#endif

	 if(mergeMode) {
		 f_out=stdout;
		 if(outfname!="stdout") {
			 f_out=fopen(outfname.chars(), "w");
			 if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
		 }
		 fprintf(f_out,"# ");
		 args.printCmdLine(f_out);
		 fprintf(f_out,"# StringTie version %s\n",VERSION);
	 }
	 else {
		 tmpfname+=".tmp";
		 f_out=fopen(tmpfname.chars(), "w");
		 if (f_out==NULL) GError("Error creating output file %s\n", tmpfname.chars());
	 }
}

//---------------
bool moreBundles() { //getter (interogation)
	bool v=true;
  v = ! NoMoreBundles;
  return v;
}

void noMoreBundles() {
	  NoMoreBundles=true;
}

void processBundle(BundleData* bundle) {
	if (verbose) {
		printTime(stderr);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->readlist.Count(), bundle->junction.Count(),
                bundle->keepguides.Count());
	#ifdef GMEMTRACE
			double vm,rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
			}
	#endif
	}
#ifdef B_DEBUG
	for (int i=0;i<bundle->keepguides.Count();++i) {
		GffObj& t=*(bundle->keepguides[i]);
		RC_TData* tdata=(RC_TData*)(t.uptr);
		fprintf(dbg_out, ">%s (t_id=%d) %s%c %d %d\n", t.getID(), tdata->t_id, t.getGSeqName(), t.strand, t.start, t.end );
		for (int fe=0;fe < tdata->t_exons.Count(); ++fe) {
			RC_Feature& exoninfo = *(tdata->t_exons[fe]);
			fprintf(dbg_out, "%d\texon\t%d\t%d\t%c\t%d\t%d\n", exoninfo.id, exoninfo.l, exoninfo.r,
					    exoninfo.strand, exoninfo.rcount, exoninfo.ucount);
			if (! (exoninfo==*(bundle->rc_data->guides_RC_exons->Get(exoninfo.id-1))))
				 GError("exoninfo with id (%d) not matching!\n", exoninfo.id);
		}
		for (int fi=0;fi < tdata->t_introns.Count(); ++fi) {
			RC_Feature& introninfo = *(tdata->t_introns[fi]);
			fprintf(dbg_out, "%d\tintron\t%d\t%d\t%c\t%d\t%d\n", introninfo.id, introninfo.l, introninfo.r,
					introninfo.strand, introninfo.rcount, introninfo.ucount);
			if (! (introninfo==*(bundle->rc_data->guides_RC_introns->Get(introninfo.id-1))))
				 GError("introninfo with id (%d) not matching!\n", introninfo.id);
		}
		//check that IDs are properly assigned
		if (tdata->t_id!=(uint)t.udata) GError("tdata->t_id(%d) not matching t.udata(%d)!\n",tdata->t_id, t.udata);
		if (tdata->t_id!=bundle->rc_data->guides_RC_tdata->Get(tdata->t_id-1)->t_id)
			 GError("tdata->t_id(%d) not matching rc_data[t_id-1]->t_id (%d)\n", tdata->t_id, bundle->rc_data->g_tdata[tdata->t_id-1]->t_id);

	}
#endif
	infer_transcripts(bundle);

	if (ballgown && bundle->rc_data) {
		rc_update_exons(*(bundle->rc_data));
	}
	if (bundle->pred.Count()>0 || ((eonly || geneabundance) && bundle->keepguides.Count()>0)) {
		if(mergeMode) GeneNo=printMergeResults(bundle, GeneNo,bundle->refseq);
		else GeneNo=printResults(bundle, GeneNo, bundle->refseq);
	}

	if (bundle->num_fragments) {
		Num_Fragments+=bundle->num_fragments;
		Frag_Len+=bundle->frag_len;
		Cov_Sum+=bundle->sum_cov;
	}

	if (verbose) {
	  /*
	  SumReads+=bundle->sumreads;
	  SumFrag+=bundle->sumfrag;
	  NumCov+=bundle->num_cov;
	  NumReads+=bundle->num_reads;
	  NumFrag+=bundle->num_frag;
	  NumFrag3+=bundle->num_fragments3;
	  SumFrag3+=bundle->sum_fragments3;
	  fprintf(stderr,"Number of fragments in bundle: %g with length %g\n",bundle->num_fragments,bundle->frag_len);
	  fprintf(stderr,"Number of fragments in bundle: %g with sum %g\n",bundle->num_fragments,bundle->frag_len);
	  */
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->pred.Count());
	#ifdef GMEMTRACE
		    double vm,rsm;
		    get_mem_usage(vm, rsm);
		    GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
		    if (rsm>maxMemRS) {
			    maxMemRS=rsm;
			    maxMemVM=vm;
			    maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		    }
	#endif
	    }
	bundle->Clear();
}





