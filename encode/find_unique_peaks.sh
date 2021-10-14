#!/bin/bash 

in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks

mkdir -p $out

for i in ${in}/*.chm13v1_peaks.bed; do
        echo $i
        base="$(basename "$i" .chm13v1_peaks.bed)"
        echo $base
	bedtools subtract -a ${in}/${base}.chm13v1_peaks.bed -b ${in}/${base}.GRCh38p13LO_peaks.bed -A > ${out}/${base}_uniquePeaks.bed
done
