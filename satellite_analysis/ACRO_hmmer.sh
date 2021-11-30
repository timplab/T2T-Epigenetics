#!/bin/bash 
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
base=ACRO
path=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO/annotate_units
# based on MSA plot pulling monomer from chr21 is best candidate
# monomer coordinates chr21:81,160-87,471
echo "chr21	81160	87471" > ${path}/${base}_Unit.coords.bed
bedtools getfasta -fi $ref -bed ${path}/${base}_Unit.coords.bed > ${path}/${base}_Unit.fasta
seq=${path}/${base}_Unit.fasta
# get  bed for the ACRO composites
annot=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CompositeRepeats.bed
grep "ACRO" $annot > ${path}/${base}_coords.bed

bedtools getfasta -fi $ref -bed ${path}/${base}_coords.bed > ${path}/${base}.fasta
fasta_formatter -i $seq a -w 80 -o ${path}/${base}_Referenceformatted.fasta #format fastas first

seqtk seq -r ${path}/${base}_Referenceformatted.fasta > ${path}/${base}_ReferenceformattedRC.fasta

muscle -in ${path}/${base}_ReferenceformattedRC.fasta -out ${path}/${base}_Referenceformatted_alnRC.fasta
muscle -in ${path}/${base}_Referenceformatted.fasta -out ${path}/${base}_Referenceformatted_aln.fasta


hmmbuild ${path}/${base}_hmm.out ${path}/${base}_Referenceformatted_aln.fasta
hmmpress  ${path}/${base}_hmm.out
nhmmscan --cpu 48 --tblout ${path}/${base}_hmmsearch.tsv ${path}/${base}_hmm.out ${path}/${base}.fasta

hmmbuild ${path}/${base}_hmmRC.out ${path}/${base}_Referenceformatted_alnRC.fasta
hmmpress  ${path}/${base}_hmmRC.out
nhmmscan --cpu 48 --tblout ${path}/${base}_hmmsearchRC.tsv ${path}/${base}_hmmRC.out ${path}/${base}.fasta


