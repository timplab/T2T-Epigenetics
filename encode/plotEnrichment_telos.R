#!/usr/bin/env Rscript

### This script takes output of encode reads aligned per repeat array and generates boxplots for log2 fold enrichment
#load libraries
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)

rm <- read_tsv(paste0(dat, "/annotations/chm13v1_polish-033121_minuschrM.bed"), col_names=F) %>%
  filter(X8 == "subtelo") %>%
  group_by(X1) %>%
  dplyr::mutate(chr_num = str_extract(X1, "[^chr]")) %>%
  mutate(rep_num=row_number()-1) %>%
  unite(id, chr_num, rep_num, sep="_")%>%
  unite(new_id, X8, id, sep="_") %>%
  rename(X1="chr", X2="start", X3="end", X4="rep_type", X5="len", X6="direction", new_id="Array", X9="perc_div") %>%
  GRanges()

ends <- read_tsv(paste0(dat, "/reference/chm13.draft_v1.0.fasta.fai"), col_names = c("chr", "end")) %>%
  mutate(start = end - 20000) %>%
  mutate(side="start")

starts=data.frame(chr=ends$chr, start=rep(0, length(ends$chr))) %>%
  mutate(start=1, end=start + 20000) %>%
  mutate(side = "end")

telos <- rbind(ends, starts) %>%
  GRanges()

telo_ovls <- FindOvls(rm,telos)

en <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/subtelomere/chm13v1_array_enrichments.tsv")

totals <- en %>%
  filter(Array == "Total") %>%
  spread(Array, Mark) %>%
  dplyr::rename("total_control"=Control, "total_treat"=Treat, "Mark"=Total)

reads <- en %>%
  filter(Array != "Total")

sum <- merge(reads, totals, by=c("Cell type", "Mark"))

stat.sum <- sum %>%
  filter(grepl("subtelo", Array)) %>%
  mutate(reg=ifelse(Array %in% telo_ovls$Array ,"telo", "non-telo")) %>%
  ungroup() %>%
  group_by(`Cell type`, reg, Mark,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) 

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich, fill=reg))+geom_boxplot(position=position_dodge(0.8), outlier.shape = NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8))+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")#+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_TAR_TELO_dotplot.pdf"),
  plot = p1,
  scale = 1,
  width = 8,
  height = 5,
)


library(ggpubr)

stat.sum <- as.data.frame(rm) %>%
  mutate(reg=ifelse(Array %in% telo_ovls$Array ,"telo", "non-telo")) %>%
  ungroup() 

p <- ggplot(stat.sum, aes(x=reg,y=perc_div))+geom_violin()+geom_dotplot(aes(fill=reg),binaxis='y', stackdir='center', dotsize=.7,position=position_dodge(.75))+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),  geom="pointrange", color="black")+labs(y="Percent Divergence", x="Region")
p + stat_compare_means(method = "kruskal.test")

ggsave(
  paste0(figs, "/", "Divergence_TAR_TELO_dotplot.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)

