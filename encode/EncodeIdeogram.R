library(BSgenome.t2t.v1.0.release)
library(tidyverse)
library(karyoploteR)
library(liftOver)
library(RColorBrewer)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

indir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks/cat"
end="unique.bed"
sample_1 <- read_tsv(paste0(indir, "/H3K27ac_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_2 <- read_tsv(paste0(indir, "/H3K27me3_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_3 <- read_tsv(paste0(indir, "/H3K36me3_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()

sample_4 <- read_tsv(paste0(indir, "/H3K4me1_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_5 <- read_tsv(paste0(indir, "/H3K4me3_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_6 <- read_tsv(paste0(indir, "/H3K9me3_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
sample_7 <- read_tsv(paste0(indir, "/CTCF_", end),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  GRanges()
annoData<-import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.combined.v4.bed")

HLA.genes <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/GOLGA.CHM13.combined.v4.bed")
novel.genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.novel.genes.v1.0.filtered.bed", col_names = c("chr", "start", "end", "name", "type")) %>%
  GRanges()

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

#chr6:31,338,352-31,642,169
detail.region <- toGRanges(data.frame("chr6", 28000000, 32000000))

for (i in seqnames(BSgenome.t2t.v1.0.release)){
  

chrom=i
colors <- brewer.pal(n = 8, name = "Set1")
win=25000
chrs=seqlevels(BSgenome.t2t.v1.0.release)
pdf(paste0(figs, "/", chrom, "_allPeaks.pdf"), width = 10, height = 8)
kp <- plotKaryotype(genome=BSgenome.t2t.v1.0.release,chromosomes=chrom,plot.type=2)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kp <- kpPlotDensity(kp, data=sample_1, data.panel=1, chr=chrom, col="#8844FF", window.size= win, r0=0, r1=0.1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=0, r1=0.1)
kp <- kpPlotDensity(kp, data=sample_3, data.panel=1,chr=chrom, col="#CC88FF", window.size= win, r0=.2, r1=.3)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.2, r1=.3)
kp <- kpPlotDensity(kp, data=sample_4,chr=chrom, data.panel=1, col="#EEAAFF", window.size= 100, r0=.4, r1=.4)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.4, r1=.5)
kp <- kpPlotDensity(kp, data=sample_5,chr=chrom,data.panel=1, col="#AA66FF", window.size= win, r0=.6, r1=.7)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.6, r1=.7)

kp <- kpPlotDensity(kp, data=sample_7,chr=chrom,data.panel=1, col="orange", window.size= win, r0=.8, r1=.9)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.8, r1=.9)

kp <- kpPlotDensity(kp, data=sample_2, chr=chrom,data.panel=1, col="red", window.size= win, r0=1, r1=1.1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2, r0=1, r1=1.1)
kp <- kpPlotDensity(kp, data=sample_6, chr=chrom,data.panel=1, col="red", window.size= win, r0=1.2, r1=1.3)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2, r0=1.2, r1=1.3)


kp <- kpPlotDensity(kp, annoData,data.panel=2, window.size= win, r0=.2, r1=.4)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density,data.panel=2, numticks = 2,r0=.2, r1=.4)
kp <- kpPlotDensity(kp, novel.genes,data.panel=2, window.size= win, r0=.4, r1=.6,col="red")
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density,data.panel=2, numticks = 2,r0=.4, r1=.6)
kpRect(kp, synteny,data.panel=2, y0=.4, y1=.6, r0=.8, r1=1,col="orange")
kpRect(kp, data= censat, y0=0, y1=1, col= colors, data.panel="ideogram", border=NA)
dev.off()

}

detail.region <- toGRanges(data.frame("chr22", 25000000, 35000000))


chrom="chr15"
colors <- brewer.pal(n = 8, name = "Set1")
win=5000
chrs=seqlevels(BSgenome.t2t.v1)
pdf(paste0(figs, "/", chrom, "GOLGA_novelPeaks.pdf"), width = 10, height = 8)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes=chrom,plot.type=2,zoom=detail.region)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kp <- kpPlotDensity(kp, data=sample_1, data.panel=1, chr=chrom, col="darkgreen", window.size= win, r0=0, r1=0.1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=0, r1=0.1)
kp <- kpPlotDensity(kp, data=sample_3, data.panel=1,chr=chrom, col="darkgreen", window.size= win, r0=.2, r1=.3)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.2, r1=.3)
kp <- kpPlotDensity(kp, data=sample_4,chr=chrom, data.panel=1, col="darkgreen", window.size= 100, r0=.4, r1=.4)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.4, r1=.5)
kp <- kpPlotDensity(kp, data=sample_5,chr=chrom,data.panel=1, col="darkgreen", window.size= win, r0=.6, r1=.7)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.6, r1=.7)

kp <- kpPlotDensity(kp, data=sample_7,chr=chrom,data.panel=1, col="orange", window.size= win, r0=.8, r1=.9)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.8, r1=.9)

kp <- kpPlotDensity(kp, data=sample_2, chr=chrom,data.panel=1, col="red", window.size= win, r0=1, r1=1.1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2, r0=1, r1=1.1)
kp <- kpPlotDensity(kp, data=sample_6, chr=chrom,data.panel=1, col="red", window.size= win, r0=1.2, r1=1.5)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2, r0=1.2, r1=1.3)

kpRect(kp, desert,data.panel=2, y0=.8, y1=1,r0=.8, r1=1, col="black")
kpRect(kp, data= censat, y0=0, y1=1, col= colors, data.panel="ideogram", border=NA)
kpPlotRegions(kp, encode_black.gr,data.panel=2, r0=.6, r1=.8,col="red")
kpPlotRegions(kp, data=HLA.genes,data.panel=2,r0=.1, r1=1)
dev.off()
