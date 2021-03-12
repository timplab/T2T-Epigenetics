#!/usr/bin/snakemake --snakefile
import pandas as pd
import os
from snakemake.utils import validate
#samples_tb = pd.read_csv(config['codedir']+"/sample_info.csv",comment="#")
repo="/home/roham/software/bismark/Bismark-0.22.2"
configfile: "snakemake_config.yml"
samples_tb = [f.split(".")[0] for f in os.listdir(config['workdir']) if f.endswith("_merge.R2.fastq.gz")]
print(str(len(samples_tb)) + " samples will be processed!")
print(samples_tb)
cores=config['threads']

##################################################
#
# alignment
#
##################################################
rule all:
  input:
   expand(config["workdir"]+"/bsseq/bismark/{sample}.R1_bismark_bt2_pe.bam", sample=samples_tb),
   expand(config['workdir']+"/bsseq/CpGcalls/{sample}", sample=samples_tb)


rule bismark_prepare_genome:
  input:
        ref = config['reference']
  output:
        out = config['workdir'] + "/Bisulfite_Genome"
  log:
   	config['workdir'] + "/bsgenomeprep.log"
  shell:"""
  	{repo}/bismark_genome_preparation --verbose {input.ref} &> {log} && touch {output.out}
     """
		

rule bismark_align_pe:
  input:
        R1=config['workdir']+"/{sample}.R1.fastq.gz",
        R2=config['workdir']+"/{sample}.R2.fastq.gz",
        bsrefdir=config['workdir']+"/Bisulfite_Genome",
        ref = config['reference']
  output:
        out = config['workdir']+"/bsseq/bismark/{sample}.R1_bismark_bt2_pe.bam"
	
  threads: config['threads']
  params:
        p=int(cores/4),
        tmpdir= config['workdir']+"/tmp/{sample}",
        outdir =  config['workdir']+"/bsseq/bismark"
  log:
        log=config['workdir']+"/bismark/{sample}.align.log"
  shell: """[ -e {params.tmpdir} ]||mkdir -p {params.tmpdir} && {repo}/bismark --bam --non_directional --bowtie2 \
		-p {params.p} --genome {input.ref} \
		-1 {input.R1} -2 {input.R2} \
		--temp_dir {params.tmpdir} \
		--output_dir {params.outdir} &> {log.log} && touch {output.out} 
	"""

rule methylation_extract:
   input:
         bam =  config['workdir']+"/bsseq/bismark/{sample}.R1_bismark_bt2_pe.bam",
         bsrefdir =  config['reference']
   output:
         outdir = config['workdir']+"/bsseq/CpGcalls/{sample}"
   shell: """[ -e {output.outdir} ]||mkdir -p {output.outdir} && {repo}/bismark_methylation_extractor -s /
   --comprehensive --merge_non_CpG \
	   -o {output.outdir} \
	   --bedGraph \
	   --remove_spaces \
	   --CX --cytosine_report \
	   --genome_folder {input.bsrefdir} \
	   {input.bam}
   """
