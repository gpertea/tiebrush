## Features to implement

*option to write a BigWig file with the coverage data
*option to write a BED junction file which can be loaded by IGV
 as an independent junction track (the _score_ field indicating depth of coverage)
 (http://software.broadinstitute.org/software/igv/splice_junctions)
 
### Note from IGV documentation about BED junction files
The splice junction view in IGV can be loaded indpendent of alignments 
by using a modified bed format, derived from the `junctions.bed"` file
produced by the TopHat program. 
 * this view is enabled by including a track line that specifies either `name=junctions` or `graphType=junctions`.
 * TopHat's "junctions.bed" file includes a track line  specifying name=junctions by default, so no action is required for these files.
