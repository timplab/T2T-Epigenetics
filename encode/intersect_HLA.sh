#!/bin/bash 

genes=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/HLA_regions.bed
tss=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.TSSonly_1kb.bed
in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_HLA_intersections
#bedtools intersect -wa -wb -a $genes -b  ${in}/*.chm13v1_peaks.bed -sorted > ${out}/novel.gene.peak.overlaps.bed

mkdir -p $out
for i in ${in}/*.chm13v1_peaks.bed; do
        echo $i
        base="$(basename "$i" .chm13v1_peaks.bed)"
        echo $base
        #bedtools subtract -a ${in}/${base}.chm13v1_peaks.bed -b ${in}/${base}.GRCh38p13LO_peaks.bed -A > ${out}/${base}_uniquePeaks.bed
	bedtools intersect -wa -wb -a $genes -b  ${i} -sorted > ${out}/NBPF_intersect_${base}.bed
done

