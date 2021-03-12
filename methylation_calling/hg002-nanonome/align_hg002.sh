#!/bin/bash

fastq=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/201126_HG002_nanonome_SRE.fastq
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/reference.fasta
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/bam
repo=/data/gmoney/Winnowmap-2.0/bin

${repo}/meryl count threads=65 k=15 output merylDB $ref
${repo}/meryl print greater-than distinct=0.9998 merylDB > ${out}/repetitive_k15.txt

${repo}/winnowmap -t 30 -W ${out}/repetitive_k15.txt -ax map-ont $ref $fastq | samtools view -Sb | \
  samtools sort -o ${out}/HG002_nanonome2-3_winnowmapk15.bam
samtools index ${out}/HG002_nanonome2-3_winnowmapk15.bam
