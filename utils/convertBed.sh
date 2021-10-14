#!/bin/bash
in="$1"
out="$2"
gff3ToGenePred -geneNameAttr=gene_name -useName $in /dev/stdout \
	     | genePredToBigGenePred /dev/stdin /dev/stdout \
	     | bedtools sort -i - > $out
