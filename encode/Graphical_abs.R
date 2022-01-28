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

mappability <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/mappability_chm13v1_200bp_score0.bed", col_names=c("chr", "start", "end", "score")) %>%
  GRanges()

gc <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome/Peaks/Cat_Peaks_BlacklistFiltered.bed", col_names = c("chr", "start", "end")) %>%
  GRanges()


annoData<-import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/CHM13.combined.v4.bed")

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

chrom=c("chr22")

colors <- brewer.pal(n = 8, name = "Set1")
win=20000
chrs=seqlevels(BSgenome.t2t.v1.0.release)
pdf(paste0(figs, "/chr22_overview.pdf"), width = 10, height = 8)
kp <- plotKaryotype(BSgenome.t2t.v1.0.release,chromosomes=chrom,plot.type=2)
kpAddBaseNumbers(kp)

kp <- kpPlotDensity(kp, novel.genes,data.panel=1, window.size= win,col="red",r0=0, r1=.2)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density,data.panel=1, numticks = 2,r0=0, r1=.2)

kp <- kpPlotDensity(kp, data=sample_1,col="#8844FF",window.size= win, r0=.3, r1=.5)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density, numticks = 2,r0=.3, r1=.5)
kpRect(kp, data= censat, y0=0, y1=1, chromosomes=chrom,col= "black", data.panel="ideogram", border=NA)
kpPlotRegions(kp, synteny,data.panel=2, r0=.1, r1=.4,col="orange")
kpPlotRegions(kp, encode_black.gr,data.panel=2, r0=.4, r1=.6, col="black")
#kpPlotRegions(kp, mappability,data.panel=2, r0=.7, r1=.8, col="red")

kp <- kpPlotDensity(kp, gc,data.panel=2, window.size= win,col="pink",r0=.8, r1=1)
max.density <- kp$latest.plot$computed.values$max.density
kpAxis(kp, ymin=0, ymax=max.density,data.panel=2, numticks = 2,r0=.8, r1=1)

dev.off()
