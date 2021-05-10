#!/usr/bin/snakemake --snakefile
import pandas as pd
import os
from snakemake.utils import validate
repo="/home/phook2/software/Bismark-0.22.2"
configfile: "snakemake_config.yml"
samples_df = pd.read_csv(config['codedir']+"/rrbs_samples.txt",names=['samples'])
samples_tb = samples_df["samples"].tolist()
print(str(len(samples_tb)) + " samples will be processed!")
print(samples_tb)
cores=config['threads']

##################################################
#
# alignment and methylation extraction
#
##################################################
rule all:
  input:
   expand(config['workdir']+"/{sample}_pass_1.fastq.gz", sample=samples_tb),
   expand(config['workdir']+"/{sample}_pass_2.fastq.gz", sample=samples_tb),
   expand(config["workdir"]+"/trimmed/{sample}_val_1.fq.gz", sample=samples_tb),
   expand(config["workdir"]+"/trimmed/{sample}_val_2.fq.gz", sample=samples_tb),
   expand(config["outdir"]+"/bsseq/{sample}/{sample}_bismark_bt2_pe.bam", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/{sample}_bismark_bt2_PE_report.txt", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/CpG_context_{sample}_bismark_bt2_pe.txt.gz", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/Non_CpG_context_{sample}_bismark_bt2_pe.txt.gz", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe_splitting_report.txt", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.M-bias.txt", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/CpG_context_{sample}_bismark_bt2_pe.txt.gz.spaces_removed.txt", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.bedGraph.gz", sample=samples_tb),
   expand(config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.bismark.cov.gz", sample=samples_tb)

rule trim_adaptors:
  input:
        R1=config['workdir']+"/{sample}_pass_1.fastq.gz",
        R2=config['workdir']+"/{sample}_pass_2.fastq.gz"
  params:
        outdir=config['workdir']+"/trimmed",
        basename="{sample}"
  log:
        log=config['workdir']+"/trimmed/{sample}_trim.log"
  output:
  	config['workdir']+"/trimmed/{sample}_val_1.fq.gz",
  	config['workdir']+"/trimmed/{sample}_val_2.fq.gz"
  shell: """[ -e {params.outdir} ]||mkdir -p {params.outdir} && trim_galore \
	  --paired \
	  --rrbs \
	  --basename {params.basename} \
	  --gzip \
	  --cores 4 \
  	  -o {params.outdir} \
	  {input.R1} {input.R2} &> {log.log}
   """

rule bismark_prepare_genome:
  input:
        ref = config['reference']
  output:
        out = config['reference'] + "/Bisulfite_Genome"
  log:
        log=config['outdir'] + "/bsgenomeprep.log"
  shell:"""
  	{repo}/bismark_genome_preparation --verbose {input.ref} &> {log.log}
     """
     
rule bismark_align_pe:
  input:
        R1_trim=config['workdir']+"/trimmed/{sample}_val_1.fq.gz",
        R2_trim=config['workdir']+"/trimmed/{sample}_val_2.fq.gz",
        ref = config['reference']
  output:
        config['outdir']+"/bsseq/{sample}/{sample}_bismark_bt2_pe.bam",
	config['outdir']+"/bsseq/{sample}/{sample}_bismark_bt2_PE_report.txt"
  threads: config['threads']
  params:
        tmpdir= config['outdir']+"/tmp/{sample}",
        outdir =  config['outdir']+"/bsseq/{sample}",
        basename = "{sample}_bismark_bt2"
  log:
        log=config['outdir']+"/bsseq/{sample}/{sample}.align.log"
  shell: """[ -e {params.tmpdir} ]||mkdir -p {params.tmpdir} && mkdir -p {params.outdir} && \
	  {repo}/bismark --bam \
	  --bowtie2 \
	  --genome {input.ref} \
	  -p 4 \
	  --basename {params.basename} \
	  -1 {input.R1_trim} -2 {input.R2_trim} \
	  --temp_dir {params.tmpdir} \
	  --output_dir {params.outdir} &> {log.log}
	"""
	
rule methylation_extract:
   input:
         bam = config['outdir']+"/bsseq/{sample}/{sample}_bismark_bt2_pe.bam",
         ref = config['reference']
   params:
         outdir = config['outdir']+"/bsseq/{sample}/CpGcalls"
   log:
         log = config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}.meth_call.log"
   output:
         config['outdir']+"/bsseq/{sample}/CpGcalls/CpG_context_{sample}_bismark_bt2_pe.txt.gz",
   	 config['outdir']+"/bsseq/{sample}/CpGcalls/Non_CpG_context_{sample}_bismark_bt2_pe.txt.gz",
   	 config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe_splitting_report.txt",
   	 config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.M-bias.txt",
    	 config['outdir']+"/bsseq/{sample}/CpGcalls/CpG_context_{sample}_bismark_bt2_pe.txt.gz.spaces_removed.txt",
   	 config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.bedGraph.gz",
         config['outdir']+"/bsseq/{sample}/CpGcalls/{sample}_bismark_bt2_pe.bismark.cov.gz"	 
   shell: """[ -e {params.outdir} ]||mkdir -p {params.outdir} && {repo}/bismark_methylation_extractor -p \
 	   --comprehensive \
           --merge_non_CpG \
	   -o {params.outdir} \
	   --bedGraph \
	   --gzip \
	   --remove_spaces \
	   --cytosine_report \
	   --genome_folder {input.ref} \
	   {input.bam} &> {log.log}
   """

