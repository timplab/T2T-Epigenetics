#!/bin/bash

array=(HG002_nanonome1_winnowmapk15 HG002_nanonome2_winnowmapk15 HG002_nanonome4-5_winnowmapk15)

for i in "${array[@]}"
do
	samtools view -h ../${i}_filtered.bam | awk 'length($10) > 20000 || $1 ~ /^@/' | samtools view | cut -f1  > ${i}_20kb_readIDs.txt
done

cat *_20kb_readIDs.txt > all_20kb_readIDs.txt

in=HG002_nanonome_CpGmethylation_allreadlengths.tsv
head -n1 $in  > head.txt
awk 'NR==FNR {a[$1]++; next}a[$5]' all_20kb_readIDs.txt $in > HG002_nanonome_CpGmethylation_20kb.tmp
cat head.txt HG002_nanonome_CpGmethylation_20kb.tmp > HG002_nanonome_CpGmethylation_20kb.tsv
rm *.tmp
~/repos/nanopolish/scripts/calculate_methylation_frequency.py -c 1.5 --split-groups HG002_nanonome_CpGmethylation_20kb.tsv > HG002_nanonome_CpGmethylationFrequency_20kb.tsv
