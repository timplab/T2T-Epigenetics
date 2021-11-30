library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(rtracklayer)
library(cowplot)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"



genes <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF.transcripts.meth.and.cutandrunMerged.bed", col_names = c("chr", "start", "end","name","n_transcripts"),delim="\t") %>%
  mutate(len=as.numeric(end)-as.numeric(start)) %>%
  mutate(ID=row_number()) %>%
  GRanges()

names=genes %>%
  as.data.frame() %>%
  dplyr::select(c(name))

indir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds/"

be2c_h3k36me3 <- read_tsv(paste0(indir,"BE2C_H3K36me3.chm13v1_peaks.bed"),col_names = c("chr", "start","end", "cell_line")) %>%
  GRanges()
be2c_h3k27me3 <- read_tsv(paste0(indir,"BE2C_H3K27me3.chm13v1_peaks.bed"),col_names = c("chr", "start","end", "cell_line"))%>%
  GRanges()

brain_microvascular_endothelial_cell_h3k36me3 <- read_tsv(paste0(indir,"brain_microvascular_endothelial_cell_H3K36me3.chm13v1_peaks.bed"),col_names = c("chr", "start","end", "cell_line"))%>%
  GRanges()

brain_microvascular_endothelial_cell_h3k27me3 <- read_tsv(paste0(indir,"brain_microvascular_endothelial_cell_H3K27me3.chm13v1_peaks.bed"),col_names = c("chr", "start","end", "cell_line"))%>%
  GRanges()

be2c_h3k36me3.dat <- FindOvls(genes,be2c_h3k36me3)%>%
  group_by(name) %>%
  mutate(npeaks=n()) %>%
  merge(names, by="name",all=TRUE) %>%
  dplyr::select(name,npeaks) %>%
  distinct() %>%
  mutate(mark="BE2C_H3K36me3")
be2c_h3k36me3.dat[is.na(be2c_h3k36me3.dat)] <- 0

be2c_h3k27me3.dat <- FindOvls(genes,be2c_h3k27me3) %>%
  group_by(name) %>%
  mutate(npeaks=n()) %>%
  merge(names, by="name",all=TRUE) %>%
  dplyr::select(name,npeaks) %>%
  distinct() %>%
  mutate(mark="BE2C_H3K27me3")
be2c_h3k27me3.dat[is.na(be2c_h3k27me3.dat)] <- 0

brain_h3k36me3.dat <- FindOvls(genes,brain_microvascular_endothelial_cell_h3k36me3,) %>%
  group_by(name)  %>%
  mutate(npeaks=n()) %>%
  merge(names, by="name",all=TRUE) %>%
  dplyr::select(name,npeaks) %>%
  distinct()%>%
  mutate(mark="brain_H3K36me3")
brain_h3k36me3.dat[is.na(brain_h3k36me3.dat)] <- 0

brain_h3k27me3.dat <- FindOvls(genes,brain_microvascular_endothelial_cell_h3k27me3) %>%
  group_by(name)  %>%
  mutate(npeaks=n()) %>%
  merge(names, by="name",all=TRUE) %>%
  dplyr::select(name,npeaks) %>%
  distinct() %>%
  mutate(mark="brain_H3K27me3")
brain_h3k27me3.dat[is.na(brain_h3k27me3.dat)] <- 0


act_all <- rbind(brain_h3k36me3.dat,be2c_h3k36me3.dat)

rep_all <- rbind(brain_h3k27me3.dat,be2c_h3k27me3.dat)

p6 <- ggplot(act_all,(aes(x=1, y=name,fill=npeaks)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "orange")+facet_wrap(~mark,ncol=4)

p7 <- ggplot(rep_all,(aes(x=1, y=name,fill=npeaks)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple")+facet_wrap(~mark,ncol=4)
plot <- plot_grid(p6,p7)

ggsave(
  paste0(figs, "/NBPF_expressionAllENCODE.pdf"),
  plot = plot,
  width = 20,
  height = 10
)

