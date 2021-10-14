#!/bin/bash
size=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chrom.sizes
gff=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.gff3


in="$1"
out="$2"
#grep "protein_coding" $gff > genes.gff 
gff3ToGenePred genes.gff -geneNameAttr=gene_name -useName /dev/stdout \
             | genePredToBigGenePred /dev/stdin /dev/stdout \
             | bedtools sort -i - > genes.bed 
awk -F "\t" genes.bed'{if($6=="+"){chr=$1;start=$2;end=$2+1;strand=$6}else{chr=$1;start=$3-1;end=$3;strand=$7};{OFS="\t"}{print chr,start,end,strand,$4}}' > CHM13.combined.v4.TSSonly.bed

bedtools slop -i CHM13.combined.v4.TSSonly.bed -g $size -b 1000 > CHM13.combined.v4.TSSonly_1kb.bed
