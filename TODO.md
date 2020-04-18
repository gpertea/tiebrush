## Features to implement

* option to write a BigWig file with the coverage data
* option to write a BED junction file which can be loaded by IGV
 as an independent junction track (the _score_ field indicating depth of coverage)
 (http://software.broadinstitute.org/software/igv/splice_junctions)
 
### Notes about BED junction files
The splice junction view in IGV can be loaded independent of alignments 
by using a modified bed format, derived from the `junctions.bed"` file
produced by the TopHat program. 
 * this view is enabled by including a track line that specifies either `name=junctions` or `graphType=junctions`.
 * TopHat's "junctions.bed" file includes a track line  specifying name=junctions by default, so no action is required for these files.
 * (from TopHat documentation): each junction line consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction
 
Example of `junctions.bed` file content:
```
track name=junctions description="TopHat junctions"
chr17_gl000205_random	36368	36643	JUNC00000003	4	+	36368	36643	255,0,0	2	21,65	0,210
chr17_gl000205_random	44994	69874	JUNC00000004	5	+	44994	69874	255,0,0	2	28,47	0,24833
chr4_gl000193_random	32293	34802	JUNC00000009	1	-	32293	34802	255,0,0	2	36,39	0,2470
chr4_gl000193_random	74344	75101	JUNC00000012	16	-	74344	75101	255,0,0	2	72,53	0,704
```
