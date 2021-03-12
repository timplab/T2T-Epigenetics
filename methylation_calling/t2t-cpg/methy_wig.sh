#!/bin/bash 
freq=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/methylation_frequency.tsv
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls
fasta=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
awk 'FNR > 1 {printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' $freq | sort -k1,1 -k2,2n > ${out}/methylation_frequency.bedgraph
bioawk -c fastx '{print $name"\t"length($seq)}' $fasta > ${out}/chrom.sizes
bedGraphToBigWig ${out}/methylation_frequency.bedgraph ${out}/chrom.sizes ${out}/methylation_frequency.bigwig
