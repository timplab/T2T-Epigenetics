library(BSgenome.t2t.v1.0.release)
library(tidyverse)
library(karyoploteR)
library(liftOver)


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

indir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks/cat"
end="unique.bed"
sample_1 <- read_tsv(paste0(indir, "/all_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

annoData<-import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.combined.v4.bed")

HLA.genes <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/GOLGA.CHM13.combined.v4.bed")

synteny <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/GRCh38_Synteny_25Kb.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
  GRanges()

novel.genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.novel.genes.v1.0.filtered.bed", col_names = c("chr", "start", "end", "name", "type")) %>%
  GRanges()

t2t.ch <- import.chain("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/bisulfite/liftover/t2t-chm13-v1.0.hg38.over.chain")
encode_black <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/ENCFF356LFX.bed", col_names = c("chr", "start", "end")) %>%
  GRanges()

encode_black.chm <- liftOver(encode_black, t2t.ch)
encode_black.gr <- unlist(encode_black.chm)



censat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
  GRanges()

desert <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDeserts.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
  GRanges()
sample_cns <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF.transcripts.meth.and.cutandrunMerged.bed",col_names = c("chr", "start", "end", "gene_id", "score")) %>%
  mutate(ID = row_number()) %>% 
  GRanges()

#chrom=c("chr1", "chr5","chr6","chr15","chr16","chr17","chr19")

colors <- brewer.pal(n = 8, name = "Set1")
win=5000
chrs=seqlevels(BSgenome.t2t.v1.0.release)
postscript(paste0(figs, "/allPeaks_chr1.eps"), width = 10, height = 8)
kp <- plotKaryotype(BSgenome.t2t.v1.0.release,chromosomes="chr1",plot.type=2)
kpAddBaseNumbers(kp)

kp <- kpPlotDensity(kp, novel.genes,data.panel="ideogram", window.size= win,col="red")
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density,data.panel="ideogram", numticks = 2)

kp <- kpPlotDensity(kp, data=sample_1,col="#8844FF",window.size= win, r0=0, r1=1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=0, r1=1)
kpRect(kp, data= censat, y0=0, y1=1, chromosomes=chrom,col= "black", data.panel="ideogram", border=NA)
kpPlotRegions(kp, synteny,data.panel=2, r0=.1, r1=.2)
kpPlotMarkers(kp, data=sample_cns,data.panel=1,r0=.1, r1=1,labels=sample_cns$gene_id,text.orientation = "horizontal",label.dist = 0.01, max.iter = 1000)
dev.off()
