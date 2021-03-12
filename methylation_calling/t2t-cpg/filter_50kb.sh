#!/bin/bash
# generate methylation frequency file for reads bigger than 50kb only 

awk 'NR==FNR {a[$1]++; next}a[$5]' 50kb_readIDs.txt methylation_calls.tsv > methylation_calls_50kb.tsv
 cat head.txt methylation_calls_50kb.tsv > methylation_calls_50kb_withhead.tsv
~/repos/nanopolish/scripts/calculate_methylation_frequency.py -c 1.5 --split-groups  methylation_calls_50kb_withhead.tsv > methylation_frequency_50kb_split.tsv
samtools view -h ont.primary.bam | awk 'length($10) > 50000 || $1 ~ /^@/' | samtools view -Sb > ont.primary50kb.bam
