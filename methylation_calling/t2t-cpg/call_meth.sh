#!/bin/bash
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta 
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls


if [ "$1" == "bed" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -q cpg -c 1.5 -i ${outdir}/methylation_calls.tsv -g $ref > ${outdir}/methylationCpG.tmp
	sort ${outdir}/methylationCpG.tmp -k1,1 -k2,2n | bgzip > ${outdir}/CpGmethylation.bed.gz
	tabix -p bed ${outdir}/CpGmethylation.bed.gz
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/CpGmethylation.bed.gz > ${outdir}/bismark.out
fi


if [ "$1" == "bam" ]; then
	dir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/bam
	base=output
	samtools view -h -b -F 272 ${outdir}/ont.primary.bam > ${outdir}/ont.primary_filtered.bam
	samtools index ${outdir}/ont.primary_filtered.bam
	~/miniconda3/bin/python /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 65 --verbose -b ${outdir}/ont.primary.bam \
		-c ${outdir}/CpGmethylation.bed.gz -f $ref |\
		samtools sort -o ${outdir}/ont.primary_meth.bam
	samtools index ${outdir}/ont.primary_meth.bam
fi
