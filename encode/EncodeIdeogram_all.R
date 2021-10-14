library(BSgenome.t2t.v1)
library(tidyverse)
library(karyoploteR)
library(liftOver)


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/sd/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks/beds"


sample_1 <- read_tsv(paste0(dat, "/HAP-1_H3K27ac.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_2 <- read_tsv(paste0(dat,"/HAP-1_H3K27me3.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_3 <- read_tsv(paste0(dat,"/HAP-1_H3K36me3.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_4 <- read_tsv(paste0(dat,"/HAP-1_H3K4me1.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_5 <- read_tsv(paste0(dat,"/HAP-1_H3K4me1.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_6 <- read_tsv(paste0(dat,"/HAP-1_H3K9me3.chm13v1_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_liftover_peaks/beds"

sample_7 <- read_tsv(paste0(dat, "/HAP-1_H3K27ac_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_8 <- read_tsv(paste0(dat,"/HAP-1_H3K27me3_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_9 <- read_tsv(paste0(dat,"/HAP-1_H3K36me3_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_10 <- read_tsv(paste0(dat,"/HAP-1_H3K4me1_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_11 <- read_tsv(paste0(dat,"/HAP-1_H3K4me1_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_12 <- read_tsv(paste0(dat,"/HAP-1_H3K9me3_GRCh38p13LO_peaks.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
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


zoom.region <- toGRanges(data.frame("chr1", 143594964, 143666019))
seqlevels(BSgenome.t2t.v1)
#pdf(paste0(figs, "/WASH_ideogram.pdf"), width = 5, height = 10)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes="chr1",plot.type=2)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kpPlotDensity(kp, data=sample_1, data.panel=1, col="#8844FF", window.size= 100000, r0=0, r1=0.1)
kpPlotDensity(kp, data=sample_2, data.panel=1, col="#AA66FF", window.size= 100000, r0=.1, r1=0.2)
kpPlotDensity(kp, data=sample_3, data.panel=1, col="#CC88FF", window.size= 100000, r0=.2, r1=.3)
kpPlotDensity(kp, data=sample_4, data.panel=1, col="#EEAAFF", window.size= 100000, r0=.3, r1=.4)
kpPlotDensity(kp, data=sample_5, data.panel=1, col="#AA66FF", window.size= 100000, r0=.4, r1=.5)
kpPlotDensity(kp, data=sample_6, data.panel=1, col="#AA66FF", window.size= 100000, r0=.5, r1=.6)

kpPlotDensity(kp, data=sample_7, data.panel=2, col="#8844FF", window.size= 100000, r0=.1, r1=.2)
kpPlotDensity(kp, data=sample_8, data.panel=2, col="#AA66FF", window.size= 100000, r0=.2, r1=0.3)
kpPlotDensity(kp, data=sample_9, data.panel=2, col="#CC88FF", window.size= 100000, r0=.3, r1=.4)
kpPlotDensity(kp, data=sample_10, data.panel=2, col="#EEAAFF", window.size= 100000, r0=.4, r1=.5)
kpPlotDensity(kp, data=sample_11, data.panel=2, col="#AA66FF", window.size= 100000, r0=.5, r1=.6)
kpPlotDensity(kp, data=sample_12, data.panel=2, col="#AA66FF", window.size= 100000, r0=.6, r1=.7)
kpPlotDensity(kp, annoData,data.panel=1, window.size= 500000, r0=.6, r1=.7)
#kpPlotRegions(kp, synteny,data.panel=1, r0=.8, r1=.9)
kpPlotRegions(kp, encode_black.gr,data.panel=1, r0=.9, r1=1)
#dev.off()





zoom.region <- toGRanges(data.frame("chr1", 143479254, 144774273))
seqlevels(BSgenome.t2t.v1)
#pdf(paste0(figs, "/WASH_ideogram.pdf"), width = 5, height = 10)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes="chr1",plot.type=2,zoom=zoom.region)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kpPlotRegions(kp, data=sample_1, data.panel=1, col="#8844FF", window.size= 100000, r0=0, r1=0.1)
kpPlotRegions(kp, data=sample_2, data.panel=1, col="#AA66FF", window.size= 100000, r0=.1, r1=0.2)
kpPlotRegions(kp, data=sample_3, data.panel=1, col="#CC88FF", window.size= 100000, r0=.2, r1=.3)
kpPlotRegions(kp, data=sample_4, data.panel=1, col="#EEAAFF", window.size= 100000, r0=.3, r1=.4)
kpPlotRegions(kp, data=sample_5, data.panel=1, col="#AA66FF", window.size= 100000, r0=.4, r1=.5)
kpPlotRegions(kp, data=sample_6, data.panel=1, col="#AA66FF", window.size= 100000, r0=.5, r1=.6)

kpPlotRegions(kp, data=sample_7, data.panel=2, col="#8844FF", window.size= 100000, r0=.1, r1=.2)
kpPlotRegions(kp, data=sample_8, data.panel=2, col="#AA66FF", window.size= 100000, r0=.2, r1=0.3)
kpPlotRegions(kp, data=sample_9, data.panel=2, col="#CC88FF", window.size= 100000, r0=.3, r1=.4)
kpPlotRegions(kp, data=sample_10, data.panel=2, col="#EEAAFF", window.size= 100000, r0=.4, r1=.5)
kpPlotRegions(kp, data=sample_11, data.panel=2, col="#AA66FF", window.size= 100000, r0=.5, r1=.6)
kpPlotRegions(kp, data=sample_12, data.panel=2, col="#AA66FF", window.size= 100000, r0=.6, r1=.7)
#kpPlotRegions(kp, annoData,data.panel=1, window.size= 500000, r0=.6, r1=.7)
#kpPlotRegions(kp, synteny,data.panel=1, r0=.8, r1=.9)
kpPlotRegions(kp, encode_black.gr,data.panel=1, r0=.9, r1=2)
#dev.off()