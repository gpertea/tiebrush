# TieBrush
Summarize and filter read alignments from multiple sequencing samples (taken as sorted BAM files).

This utility aims to merge/collapse "duplicate" read alignments (same location with the same CIGAR string), across multiple sequencing samples (multiple input BAM files), adding custom SAM tags in order to keep track of the "alignment multiplicity" count (how many times the same alignment is seen across all input data) and "sample count" (how many samples show that same alignment).

The initial goal is to generate this composite BAM file which multiplexes read alignments from many sequencing samples, painting a comprehensive "background" picture of read alignments with their counts across many samples.

(some alignment filtering can also be implemented if needed).

# SAM tags implemented
* __YC__:i:nn : nn = how many alignments were collapsed into an alignment record (duplicity count)
* __YX__:i:nn : nn = how many samples have this alignment (sample count)

If either of these tags are missing (i.e. GBamRecord::__tag_int__() call returns 0) then the alignment is unique (when YC is 0) or only one sample has it (if YX is 0). The actual count in these cases is obviously 1. 

