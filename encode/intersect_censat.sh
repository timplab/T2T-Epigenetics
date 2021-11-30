#!/bin/bash 

censat=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations/t2t-chm13.v1.0.cenSat_regions_sorted.bed
in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/H3K9me3_intersections
mkdir -p $out
for i in ${in}/*H3K9me3.chm13v1_peaks.bed; do
        echo $i
        base="$(basename "$i" .chm13v1_peaks.bed)"
        echo $base
        #bedtools subtract -a ${in}/${base}.chm13v1_peaks.bed -b ${in}/${base}.GRCh38p13LO_peaks.bed -A > ${out}/${base}_uniquePeaks.bed
        bedtools intersect -a $i -b  $censat -sorted > ${out}/censat_intersect_${base}.bed
	bedtools intersect -v -a $i -b  $censat -sorted > ${out}/censat_nonintersect_${base}.bed
done
