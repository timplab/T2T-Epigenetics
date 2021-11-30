library(BSgenome.t2t.v1)
library(tidyverse)
library(karyoploteR)
library(liftOver)
library(RColorBrewer)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks/beds"

tail="_H3K9me3.chm13v1_peaks.bed"

sample_1 <- read_tsv(paste0(dat, "/BE2C", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_2 <- read_tsv(paste0(dat, "/Caco-2", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_3 <- read_tsv(paste0(dat, "/HAP-1", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_4 <- read_tsv(paste0(dat, "/MG63", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_5 <- read_tsv(paste0(dat, "/SJCRH30", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_6 <- read_tsv(paste0(dat, "/SJSA1", tail),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

annoData<-import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.bed")

synteny <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/GRCh38_Synteny_25Kb.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
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

colors <- brewer.pal(n = 8, name = "Set1")
win=10000
seqlevels(BSgenome.t2t.v1)
pdf(paste0(figs, "/ACRO_H3K9me3.pdf"), width = 5, height = 20)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes=c("chr13","chr14","chr15", "chr21", "chr22"),plot.type=2)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kpPlotDensity(kp, data=sample_1, data.panel=1, col="#8844FF", window.size= win, r0=0, r1=0.1)
kpPlotDensity(kp, data=sample_2, data.panel=1, col="#AA66FF", window.size= win, r0=.1, r1=0.2)
kpPlotDensity(kp, data=sample_3, data.panel=1, col="#CC88FF", window.size= win, r0=.2, r1=.3)
kpPlotDensity(kp, data=sample_4, data.panel=1, col="#EEAAFF", window.size= win, r0=.3, r1=.4)
kpPlotDensity(kp, data=sample_5, data.panel=1, col="#AA66FF", window.size= win, r0=.4, r1=.5)
kpPlotDensity(kp, data=sample_6, data.panel=1, col="#AA66FF", window.size= win, r0=.5, r1=.6)
#kpPlotDensity(kp, annoData,data.panel=1, window.size= win, r0=.6, r1=.7)
kpPlotRegions(kp, synteny,data.panel=1, r0=.8, r1=.9)
#kpPlotRegions(kp, encode_black.gr,data.panel=1, r0=.9, r1=1)
#kpRect(kp, desert,data.panel=1, y0=1, y1=1.5,r0=1, r1=1.1, col="black")
kpRect(kp, data= censat, y0=0, y1=1, col= colors, data.panel="ideogram", border=NA)
dev.off()


