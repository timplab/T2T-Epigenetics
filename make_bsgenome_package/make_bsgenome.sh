#!/bin/bash 
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta

out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/seqnames
# genome needs to be split by chr 

for i in {1..22};
do
	samtools faidx ${ref} chr${i} > ${out}/chr${i}.fa
done

samtools faidx ${ref} chrX > ${out}/chrX.fa
samtools faidx ${ref} chrM > ${out}/chrM.fa
