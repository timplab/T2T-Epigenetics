#!/bin/bash 
path=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/hsat
seq=${path}/hsat_1.77kb.fasta
base=hsat
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
# first cat together DXZ4 repeat from chm13 and the DXZ4 clone sequence from GenBank 
# align them as a msa with muscle then convert into a stockholm format for hmmer
echo "chr16	49163529	49239753" > ${path}/${base}.coords.bed		
#bedtools getfasta -fi $ref -bed ${path}/${base}.coords.bed > ${path}/${base}.fasta

fasta_formatter -i $seq a -w 80 -o ${path}/${base}_Referenceformatted.fasta format fastas first
#clustalo --in=$seq --out=${path}/msa_DXZ4_clustalo.sto -t DNA --force --outfmt=st --wrap=80 -v --threads 48

# also get reverse complement of reference se
seqtk seq -r ${path}/${base}_Referenceformatted.fasta > ${path}/${base}_ReferenceformattedRC.fasta

muscle -in ${path}/${base}_ReferenceformattedRC.fasta -out ${path}/${base}_Referenceformatted_alnRC.fasta
muscle -in ${path}/${base}_Referenceformatted.fasta -out ${path}/${base}_Referenceformatted_aln.fasta


hmmbuild ${path}/${base}_hmm.out ${path}/${base}_Referenceformatted_aln.fasta
hmmpress  ${path}/${base}_hmm.out
nhmmscan --cpu 48 --tblout ${path}/${base}_hmmsearch.tsv ${path}/${base}_hmm.out ${path}/${base}.fasta

hmmbuild ${path}/${base}_hmmRC.out ${path}/${base}_Referenceformatted_alnRC.fasta
hmmpress  ${path}/${base}_hmmRC.out
nhmmscan --cpu 48 --tblout ${path}/${base}_hmmsearchRC.tsv ${path}/${base}_hmmRC.out ${path}/${base}.fasta
# load output into R to parse with hmmer parser package
# now put together into multi fasta and do a msa
#bedtools getfasta -fi ../ref/t2t-chm13.20200727.fasta -bed DXZ4_nhmmer_parsed.tsv -split -name > DXZ4_monomers.fasta
