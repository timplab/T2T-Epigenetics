#!/bin/bash 
in=methylation_calls_50kb_withhead.tsv

awk '{A[$6]++}END{for(i in A)print i,A[i]}' $in
