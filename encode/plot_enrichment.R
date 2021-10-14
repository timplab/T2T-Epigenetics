#!/usr/bin/env Rscript

### This script takes output of encode reads aligned per repeat array and generates boxplots for log2 fold enrichment
#load libraries
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)

# load data
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
en <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/subtelomere/chm13v1_array_enrichments.tsv")

totals <- en %>%
  filter(Array == "Total") %>%
  spread(Array, Mark) %>%
  dplyr::rename("total_control"=Control, "total_treat"=Treat, "Mark"=Total)

reads <- en %>%
  filter(Array != "Total")

sum <- merge(reads, totals, by=c("Cell type", "Mark"))

SAT = c("HSAT1", "HSAT2", "HSAT3","HOR", "CT", "BSAT", "GSAT", "HSAT4", "DHOR", "MON", "subtelo-TAR")

sum <- sum %>%
  mutate(name = ifelse(grepl("hsat1", Array), "HSAT1", Array)) %>%
  mutate(name = ifelse(grepl("hsat2", Array), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("hsat3", Array), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("hsat4", Array), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("hsat5", Array), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("hsat6", Array), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", Array), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", Array), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("dhor", Array), "DHOR", name)) %>%
  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", Array), "MON", name)) %>%
  mutate(name = ifelse(grepl("GSAT", Array), "GSAT", name)) %>%
 # mutate(name = ifelse(grepl("tar", Array), "centro-TAR", name)) %>%
  mutate(name = ifelse(grepl("subtelo", Array), "subtelo-TAR", name)) %>%
  rename(Array="rep_type", name="Array") 

stat.sum <- sum %>%
  filter(Array %in% SAT) %>%
  group_by(`Cell type`, Mark,Array,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) 


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich, color=Array))+ geom_violin()+
  geom_dotplot(binaxis='y', stackdir='center')+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+facet_wrap(~Array)+scale_color_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_perCell_dotplot.pdf"),
  plot = p1,
  scale = 1,
  width = 15,
  height = 10,
)


stat.sum <- sum %>%
  group_by(`Cell type`, Mark,Array,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) %>%
  filter(Mark %in% c("H3K9me3")) %>%
  filter(Array %in% SAT)

p1 <- ggplot(stat.sum, aes(x=Array, y=enrich, fill=Array))+geom_boxplot()+ 
  geom_dotplot(binaxis='y', stackdir='center')+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+facet_wrap(~Mark, ncol=1)+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_h3k9me3_satboxplot.pdf"),
  plot = p1,
  scale = 1,
  width = 8,
  height = 5,
)

p3 <- ggplot(stat.sum, aes(x=`Cell type`, y=enrich))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = .8)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")

ggsave(
  paste0(figs, "/", "ENCODE_perCellType_dotplot.pdf"),
  plot = p3,
  scale = 1,
  width = 8,
  height = 5,
)


stat.sum <- sum %>%
  group_by(`Cell type`, Mark,Array,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) %>%
  #filter(Mark %in% c("H3K9me3")) %>%
  filter(Array %in% c("centro-TAR", "subtelo-TAR"))

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich, fill=Array))+geom_boxplot(position=position_dodge(0.8), outlier.shape = NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8))+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_TAR_dotplot.pdf"),
  plot = p1,
  scale = 1,
  width = 8,
  height = 5,
)

stat.sum <- sum %>%
  filter(Array %in% c("centro-TAR", "subtelo-TAR")) %>%
  group_by(`Cell type`, rep_type, Mark,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat)))

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich, fill=Array))+geom_boxplot(position=position_dodge(0.8), outlier.shape = NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8))+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+scale_fill_manual(values=censatColors)


