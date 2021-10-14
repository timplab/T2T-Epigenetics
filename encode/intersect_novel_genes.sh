#!/bin/bash 

genes=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.novel.genes.v1.0.filtered.bed
tss=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.TSSonly_1kb.bed
in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all
#bedtools intersect -wa -wb -a $genes -b  ${in}/*.chm13v1_peaks.bed -sorted > ${out}/novel.gene.peak.overlaps.bed
bedtools intersect -wa -wb -a $tss -b  ${in}/*H3K27ac_uniquePeaks.bed -sorted > ${out}/novel_H3K27ac.peaks.gene.overlaps.bed
