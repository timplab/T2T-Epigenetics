#!/usr/bin/env Rscript

### This script takes output of encode reads aligned per repeat array and generates boxplots for log2 fold enrichment
#load libraries
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)

# load data
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
en <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/210308_alignments/bt2.chm13v1_array_enrichmentsv2.tsv")

totals <- en %>%
  filter(Array == "Total") %>%
  spread(Array, Mark) %>%
  dplyr::rename("total_control"=Control, "total_treat"=Treat, "Mark"=Total)

reads <- en %>%
  filter(Array != "Total")

sum <- merge(reads, totals, by=c("Cell type", "Mark"))

SAT = c("HSAT1", "HSAT2", "HSAT3","HOR", "CT")

stat.sum <- sum %>%
  group_by(`Cell type`, Mark,Array,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) %>%
  filter(Mark %in% c("H3K9me3")) %>%
  filter(Array %in% SAT)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich))+ 
  geom_dotplot(binaxis='y', stackdir='center')+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+facet_wrap(~Array)

ggsave(
  paste0(figs, "/", "ENCODE_perCell_dotplot.pdf"),
  plot = p1,
  scale = 1,
  width = 8,
  height = 5,
)

p2 <- ggplot(stat.sum, aes(x=Array, y=enrich, fill=Array))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = 1.5)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_perSat_dotplot.pdf"),
  plot = p2,
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
