#!/bin/bash 
indir=/mithril/Data/T2T_Data/hg002_aws/hg002/meth_calls
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled

for file in ${indir}/*_methylation.tsv
do
	base="$(basename -- $file ".tsv")"
	echo $base
       	grep "chrX" $file > ${outdir}/${base}_chrX.tsv &
done
cat ${outdir}/head.txt ${outdir}/GM24385*.tsv  > ${outdir}/HG002_chrX_CpG_methylation.tsv
~/repos/nanopolish/scripts/calculate_methylation_frequency.py -c 1.5 --split-groups ${outdir}/HG002_chrX_CpG_methylation.tsv > ${outdir}/HG002_HPRC_pooledCpGmethylationFrequency_HG002chrX.tsv
