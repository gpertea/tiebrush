# TieBrush

Summarize and filter read alignments from multiple sequencing samples
(taken as sorted BAM files).

This utility aims to merge/collapse "duplicate" read alignments (same
location with the same CIGAR string), across multiple sequencing samples
(multiple input BAM files), adding custom SAM tags in order to keep
track of the "alignment multiplicity" count (how many times the same
alignment is seen across all input data) and "sample count" (how many
samples show that same alignment).

The initial goal is to generate this composite BAM file which
multiplexes read alignments from many sequencing samples, painting a
comprehensive "background" picture of read alignments with their counts
across many samples.

(some alignment filtering can also be implemented if needed).

# SAM tags implemented
* __YC__:i:N stores the number of alignments that were merged into this alignment record (multiplicity count)
* __YX__:i:N stores the number of samples that have this alignment (sample count)

If either of these tags are missing (i.e. GBamRecord::__tag_int__() call returns 0) then the alignment is unique (when YC is 0) or only one sample has it (if YX is 0). The actual count in these cases is obviously 1. 

In an attempt to reflect a measure of "mosaicity" for read alignments in a sample, the following tag is also added :

* __YD__:i:N keeps track of the maximum number of contiguous bases preceding the start of the read alignment in the samples(s) that it belongs to. In other words, if the current alignment is part of an exon-overlapping bundle (strand specific!), this value holds the maximum distance from the beginning of the bundle to the start of this alignment, across all samples having this alignment. If the alignment is not in a bundle (i.e. it is preceded by a uncovered region as it is not overlapped by any another alignment with a lower start position), in all the individual samples where that alignment is present, then the YD value is 0 and the tag is omitted from the output file produced by TieBrush. That means that all the alignments lacking a YD tag in the TieBrush output start at the very beginning of an exon-overlapping bundle (i.e. are not overlapped by a preceding alignment with a lower start coordinate).


# TieCov

The tiecov utility can take the output file produced by TieBrush and can generate the following auxiliary base/junction coverage files:
   * a BedGraph file with the coverage data (see http://genome.ucsc.edu/goldenPath/help/bedgraph.html); this file can be converted to BigWig (using bedGraphToBigWig) or to TDF format (using igvtools) in order to be loaded in IGV as an additional coverage track
   * a junction BED file which can be loaded directly in IGV as an additional junction track (http://software.broadinstitute.org/software/igv/splice_junctions)
