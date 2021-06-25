# Simulated Example Data

This directory contains a small simulated dataset for testing and getting 
familiar with TieBrush.

## Running TieBrush and TieCov

Output of the commands below is provided along with the simulated data and includes the following files:

1. test/t1/t1.bam and test/t2/t2.bam - collapsed representations of simulated tissue #1 and #2 computed with TieBrush
2. test/t1/t1.coverage.bedgraph and test/t2/t2.coverage.bedgraph - Coverage computed in BED/BedGraph format
3. test/t1/t1.sample.bedgraph and test/t2/t2.sample.bedgraph - Approximated number of samples computed for each mapped position in BED/BedGraph format
4. test/t1/t1.junctions.bed and test/t2/t2.junctions.bed - Splice Junction coordinates and coverage computed in BED format

### Collapsing reads and building representative alignment for each tissue with TieBrush

```tiebrush -o t1/t1.bam t1/t1s0.bam t1/t1s1.bam t1/t1s2.bam t1/t1s3.bam t1/t1s4.bam t1/t1s5.bam t1/t1s6.bam t1/t1s7.bam t1/t1s8.bam t1/t1s9.bam```

```tiebrush -o t2/t2.bam t2/t2s0.bam t2/t2s1.bam t2/t2s2.bam t2/t2s3.bam t2/t2s4.bam t2/t2s5.bam t2/t2s6.bam t2/t2s7.bam t2/t2s8.bam t2/t2s9.bam```

### Extracting BED-formatted summaries for each tissue with TieCov

```tiecov -s t1/t1.sample -c t1/t1.coverage -j t1/t1.junctions t1/t1.bam```

```tiecov -s t2/t2.sample -c t2/t2.coverage -j t2/t2.junctions t2/t2.bam```

### Visualizing results in IGV
Upon completion, tiecov will have generated several summary files for each tissue.
Along with the collapsed alignment, these summaries can be viewed and scrutinized using
genome browsers such as [IGV](http://software.broadinstitute.org/software/igv/)
that support BAM (for collapsed alignment), BED (for the junctions track)
and BEDgraph (for the sample and coverage tracks) formats.

## Simulation Protocol
1. 10 samples were randomly selected from brain tissue and another 10 samples were randomly selected 
from the heart tissue of the GTEx project [2] (Table 1). 
2. Samples were aligned with HISAT2 [3] against the GRCh.38(patch 12) [4].
3. Transcriptomes for each sample were individually assembled and quantified using StringTie2 [5].
4. Using GffCompare [6] samples were individually compared against the RefSeq annotation [7].
5. Transcript structures (GTF) and coverage data of isoforms overlapping the NEFL and SLC25A3 genes 
   were extracted for each sample.
7. Extracted transcriptome (FASTA sequences obtained with gffread [6]) and coverage data 
   were used to simulate 101bp reads using polyester [8].
8. Fasta files produced by polyester were directly converted into SAM format using sim2sam utility [9].

| Brain Samples | Heart Samples |
|---------------|:--------------|
|SRR1368222     | SRR1345221    |
|SRR1370173  	| SRR1441464    |
|SRR1374221 	| SRR1432366    |
|SRR1363810     | SRR661229     |
|SRR1467882     | SRR1337388    |
|SRR1471773  	| SRR1458484    |
|SRR1473590   	| SRR816358     |
|SRR814989   	| SRR1352213    |
|SRR817658   	| SRR1435584    |
|SRR818146   	| SRR1345736    |

**Table 1.** *GTEx samples used to model parameters for the simulation*

