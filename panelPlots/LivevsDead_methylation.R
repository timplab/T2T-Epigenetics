# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/ilee_plot_utils.R")
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/methylation_R_utils.R")
library(entropy)


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

hor = read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.HOR_annotations.bed"), col_names = c("chr", "start", "end", "name", "blank1", "blank2", "score", "strand", "nm1", "nm2")) %>%
  dplyr::filter(grepl("hor", name))%>%
  mutate(status = ifelse(grepl("L", name), "live", "dead")) %>%
  dplyr::select(chr,start,end,name,status) %>%
  mutate(width=end-start) %>%
  group_by(name) %>%
  mutate(len=sum(width)) %>%
  GRanges()

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

ovls.meth <- FindOvls(hor,chm13_meth) %>%
  dplyr::filter(!grepl("dhor", name)) %>%
  group_by(name,status,seqnames,len) %>%
  summarise(methylated_frequency=mean(methylated_frequency))

ggplot(ovls.meth, aes(x=len/1e6,y=methylated_frequency, color=status))+geom_point(size=4,alpha=.5)+geom_text(aes(label=seqnames),hjust=0, vjust=0)

ggsave(
  paste0(figs, "/LivevsDeadHOR.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)

binned_meth <-read_tsv(paste0(dat, "/chm13_final_beds/HG002_autosomes_AllCen_10kbBinned_ALL.bed")) %>%
  GRanges

cen.meth <- FindOvls(hor,chm13_meth) %>%
  group_by(name) %>%
  mutate(x = ceiling(row_number()/10)) %>%
  ungroup() %>%
  group_by(x,name,seqnames,status) %>%
  summarise(var=var(methylated_frequency), mean_meth=mean(methylated_frequency)) %>%
  group_by(name,seqnames,status) %>%
  summarise(mean_meth=mean(mean_meth), mean_ent=mean(var))

cen.meth$mean_ent[is.na(cen.meth$mean_ent)] <- 0
ggplot(cen.meth, aes(x=mean_ent,y=mean_meth, color=status))+geom_point()+geom_text(aes(label=seqnames),hjust=0, vjust=0)

mm.chm13 <- chm13_meth %>%
  as.data.frame() %>%
  group_by(seqnames) %>%
  mutate(x = ceiling(row_number()/50)) %>%
  group_by(seqnames,x) %>%
  summarise(mean_meth=mean(methylated_frequency), ent=var(methylated_frequency),start=min(start), end=max(end)) %>%
  GRanges()

cen.var <- FindOvls(hor,binned_meth) %>%
  mutate(len=end-start) %>%
  group_by(name,len,status,seqnames) %>%
  summarise(mean_meth=mean(freq), var=var(freq)/len) %>%
  distinct()

ggplot(cen.var, aes(x=mean_meth,y=var, color=status))+geom_point()+geom_text(aes(label=seqnames),hjust=0, vjust=0)
