#!/bin/bash 
repo=/home/isac/Software/nanopolish-cpggpc
bam=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/bam/whole_genome/HG002_nanonome1_winnowmapk15TELOCLUSTERS.bam 
path=/pym/Data/Nanopore/projects/hg002_nanonome
fastq=${path}/201119_HG002_nanonome_SREXL.fastq
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/reference/reference.fasta
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls
basename=HG002_nanonome1_winnowmapk15TELOCLUSTERS

samtools view -h -b -F 256 ${bam}  > ${outdir}/${basename}_filtered.bam

bam=${outdir}/${basename}_filtered.bam
samtools index $bam
${repo}/nanopolish call-methylation --bam $bam --reads $fastq --genome $ref --methylation gpc -t 90 --progress > ${outdir}/${basename}_GpCmethylation.tsv
${repo}/nanopolish call-methylation -b ${bam} -r $fastq -g $ref -q cpg -t 90 --progress > ${outdir}/${basename}_CpGmethylation.tsv

python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -q cpg -c 1.5 --nome -i ${outdir}/${basename}_CpGmethylation.tsv -g $ref > ${outdir}/${basename}_methylationCpG.tmp
python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -q gpc -c 1.0 --nome -i ${outdir}/${basename}_GpCmethylation.tsv -g $ref > ${outdir}/${basename}_methylationGpC.tmp
sort ${outdir}/${basename}_methylationCpG.tmp -k1,1 -k2,2n | bgzip > ${outdir}/${basename}_CpGmethylation.bed.gz
tabix -p bed ${outdir}/${basename}_CpGmethylation.bed.gz
sort ${outdir}/${basename}_methylationGpC.tmp -k1,1 -k2,2n | bgzip > ${outdir}/${basename}_GpCmethylation.bed.gz
tabix -p bed ${outdir}/${basename}_GpCmethylation.bed.gz
