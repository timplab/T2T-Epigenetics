#!/bin/bash 
freq=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/compare_hg002/methylation_frequency_HG002_split.tsv
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/compare_hg002
fasta=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/chm13_HG002X_HG38Y.fasta
awk 'FNR > 1 {printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3+1,$7}' $freq | sort -k1,1 -k2,2n > ${out}/methylation_frequency.bedgraph
bioawk -c fastx '{print $name"\t"length($seq)}' $fasta > ${out}/chrom.sizes
bedGraphToBigWig ${out}/methylation_frequency.bedgraph ${out}/chrom.sizes ${out}/HG002_methylation_frequency.bigwig
