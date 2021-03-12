#!/bin/bash 

in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/methylation_calls_50kb_withhead.tsv
censat=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations/t2t-chm13.v1.0.cenSat_annotationFormatted.bed
censatv2=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$3,$4,$2,$5,$6}' $in | bedtools intersect -a stdin -b $censatv2 -wb
