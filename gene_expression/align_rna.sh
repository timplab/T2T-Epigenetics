#!/bin/bash 


if [ "$1" == "chm13" ]; then
	in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/reads
	out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/bam
	base=CHM13
	asm=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
	hisat2-build $asm /kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0
	ind=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0	

	paired_files=$(ls "${in}/CHM13"* | grep R1)
	
	for i in $paired_files; do
		f=$(basename $i _R1_001_cutadapt-m100.fastq.gz)
		end=_001_cutadapt-m100.fastq.gz
		echo $f
		echo ${in}/${f}_R1${end}
		hisat2 -x $ind -1 ${in}/${f}_R1${end} -2 ${in}/${f}_R2${end} -p 72 | samtools view -Sb | samtools sort -o ${out}/${f}.bam &> ${out}/${f}_alignment.log
	done
fi
if [ "$1" == "hg002" ]; then
	in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/reads
	out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/bam
	base=CHM13
	asm=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/reference.fasta
	hisat2-build $asm /kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/reference
	ind=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/reference

	paired_files=$(ls "${in}/HG002"* | grep R1)
	
	for i in $paired_files; do
		f=$(basename $i _R1_cutadapt-m100.fastq.gz)
		end=_cutadapt-m100.fastq.gz
		echo $f
		echo ${in}/${f}_R1${end}
		hisat2 -x $ind -1 ${in}/${f}_R1${end} -2 ${in}/${f}_R2${end} -p 72 | samtools view -Sb | samtools sort -o ${out}/${f}.bam &> ${out}/${f}_alignment.log
	done
fi

