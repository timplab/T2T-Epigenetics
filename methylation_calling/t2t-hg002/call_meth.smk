"""
This Snakemake pipeline is intended to absorb HG002 pangenome sequencing data and methylation call it on AWS.  Based on Ariel's original snakemake
    * Aligns to the reference genome with winnowmap
    * Indexes methylation with f5c
    * Calls methylation with Nanopolish
    * Formats methylation calls with nanopore-methylation-utilities
    fast5 directory and fastq file must have same basename, summary file not required 
"""
configfile: "config.yaml"
###-------input paths/files -------###

import glob
import re

##WT - note check this to figure out parallelization properly
##cores=config["cores"]

##Software config
nanopolish = config["nanopolish"]
winnowmap = config["winnowmap"]
util = config["nanopore-methylation-utilities"]

##Pass in directory we are working in --config base=dir - default is /data
base=config["base"]
outdir=config["outdir"]
##Pass in which references using --config ref=chm13 or ref=hg002 - default is chm13.
reflist = config["reflist"]
   
###--------------###

##WT - again way more complicated than necessary. I noted that there is a pattern of GM24385_{0..22} I can use, broken only by {NB}
##Find the GM24385 and then get the "next" thing after that in split?
filelist=glob.glob(base+"/*GM24385*gz")
print(filelist)
idlist=[]
for filey in filelist:
    nameparts=re.split("[_\/]", filey)
    idlist.append("GM24385_"+nameparts[nameparts.index("GM24385")+1])

print(idlist)

# ###------- Pipeline Rules -------#####

workdir: base

rule all:
    input:
        expand( outdir + "/{ref}/meth_calls/{id}_meth.bam", id=idlist, ref=reflist)
    shell:
        "echo {base}"

##Make reference index        
rule ref_index:
    output:
        "{ref}.k15.txt"
    threads: workflow.cores/2-2
    run:
        reffa=reflist[wildcards.ref]
        shell("{winnowmap}/meryl count threads={threads} k=15 output {wildcards.ref}.merylDB {reffa}")
        shell("{winnowmap}/meryl print greater-than distinct=0.9998 {wildcards.ref}.merylDB > {output}")

##Align with winnowmap
rule align:
    input:
        refidx=rules.ref_index.output
    output:
        bam = outdir + "/{ref}/bam/{id}.bam"
    threads: workflow.cores/2-2
    message: """Aligning to reference with winnowmap"""
    run:
        fq=glob.glob(base+"/*"+wildcards.id+"*fastq.gz")[0]
        reffa=reflist[wildcards.ref]
        shell("{winnowmap}/winnowmap -t {threads} -W {input.refidx} -ax map-ont {reffa} {fq} |"+
              "samtools view -b -u -F 256 | samtools sort -o {output.bam}")
        shell("samtools index {output.bam}")

#rule np_index:
##In principle I need this, but I already indexed so I'm not going to bother for now.
        
rule nanopolish:
    input:
        bam = rules.align.output.bam
    output:
        meth=outdir + "/{ref}/meth_calls/{id}_CpG_methylation.tsv"
    message: """calling methylation with nanopolish"""
    threads: workflow.cores/2-2
    run:
        reffa=reflist[wildcards.ref]
        fq=glob.glob(base+"/*"+wildcards.id+"*fastq.gz")[0]
    	shell("{nanopolish}/nanopolish call-methylation -b {input.bam} -r {fq} -g {reffa} -q cpg"+
              " -t {threads} --progress > {output.meth}")

#isac methylbed code
rule methylbed:
    input:
        tsv = outdir + "/{ref}/meth_calls/{id}_CpG_methylation.tsv"
    output:
        methbed=outdir + "/{ref}/meth_calls/{id}_CpG_methylation.bed.gz"
    message: """format methyl bed"""
    threads: 1
    run:
        reffa=reflist[wildcards.ref]
        temp=outdir+"/"+wildcards.ref+"/meth_calls/"+wildcards.id+".tmp"
        shell("python3 {util}/mtsv2bedGraph.py -q cpg -c 1.5 -g {reffa} -i {input.tsv} > {temp}")
        shell("sort {temp} -k1,1 -k2,2n | bgzip > {output.methbed}")
        shell("tabix -p bed {output.methbed}")
        shell("rm {temp}")

#isac methylbam code
rule methylbam:
    input:
        bam = outdir + "/{ref}/bam/{id}.bam",
        tsv = rules.methylbed.output.methbed
    output:
        outdir + "/{ref}/meth_calls/{id}_meth.bam"
    message: """convert bam for methylation"""
    threads: workflow.cores/2-2
    run:
        reffa=reflist[wildcards.ref]
        shell("python3 {util}/convert_bam_for_methylation.py -t {threads} --windowsize 1000000 --verbose -b {input.bam}"+
              " -c {input.tsv} -f {reffa} | samtools sort -o {output}")
        shell("samtools index {output}")

# ###---###
