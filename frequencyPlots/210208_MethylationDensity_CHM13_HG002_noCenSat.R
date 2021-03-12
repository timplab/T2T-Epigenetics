#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(bsseq)
library(Biostrings)
library(ggplot2)
library(png)
library(cowplot)
options(scipen=999)
library(zoo)
library(BSgenome.t2t.v1.0.release)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

# set data output path

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/violin"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  mutate(auto = case_when(chromosome == "chrX" ~ "chrX", 
                          TRUE ~ "autosome")) %>%
  GRanges()


HG002_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/pooled/HG002_nanonome_CpGmethylationFrequency.tsv")) %>%
  filter(chromosome  != "chrY") %>%
  mutate(auto = case_when(chromosome == "chrX" ~ "chrX", 
                          TRUE ~ "autosome")) %>%
  GRanges()

# find overlaps between censat regions and methylation calls
cen <- read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_regions.bed")) %>%
  dplyr::rename("chr"=1, "start"=2, "end"=3, "name"=4) %>%
  GRanges()

keepi <- findOverlaps(chm13_meth,cen)
chm13_meth.gr <- as.data.frame(chm13_meth[-queryHits(keepi)])

keepi <- findOverlaps(HG002_meth,cen)
HG002_meth.gr <- as.data.frame(HG002_meth[-queryHits(keepi)])



violin <- ggplot(chm13_meth.gr, aes(x = auto, y = methylated_frequency, fill = auto))+geom_violin(adjust=2)+theme_classic(base_size = 20)+labs(x = "Repeat", y = "CpG Density")+geom_boxplot(outlier.shape = NA, width=.1)


ggsave(
  paste0(figs, "/chm13WG_methylationFreqNoCenSat.pdf"),
  plot = violin,
  scale = 1,
  width = 5,
  height = 5,
)

violin2 <- ggplot(HG002_meth.gr, aes(x = auto, y = methylated_frequency, fill = auto))+geom_violin(adjust=2)+theme_classic(base_size = 20)+labs(x = "Repeat", y = "CpG Density")+geom_boxplot(outlier.shape = NA, width=.1)

ggsave(
  paste0(figs, "/HG002WG_methylationFreqNoCenSat.pdf"),
  plot = violin2,
  scale = 1,
  width = 5,
  height = 5,
)
