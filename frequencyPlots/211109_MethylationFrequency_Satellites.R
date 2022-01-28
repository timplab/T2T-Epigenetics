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
library(ggridges)
library(BSgenome.t2t.v1.0.release)
options(knitr.duplicate.label = 'allow')
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
library("ggsci")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")

hg002_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv") %>%
  GRanges()

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

hg002.flag <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Blacklist.bed", col_names=c("chr", "start", "end")) %>%
  GRanges()


# set data output path

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/violin"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load and parse censat files
list=c("HSAT1", "HSAT2", "HSAT3", "HSAT4", "BSAT", 'GSAT', "ACRO","SST1_Composite", "MON","HOR", "DHOR")

censat.widths <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
  mutate(name = ifelse(grepl("hsat1", name), "HSAT1", name)) %>%
  mutate(name = ifelse(grepl("hsat2", name), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("hsat3", name), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("hsat4", name), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("hsat5", name), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("hsat6", name), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", name), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", name), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("dhor", name), "DHOR", name)) %>%
  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", name), "MON", name)) %>%
  mutate(name = ifelse(grepl("GSAT", name), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("ACRO", name), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("SST1_Composite", name), "SST1_Composite", name)) %>%
  mutate(width=chromEnd-chromStart) %>%
  filter(name %in% list) %>%
  dplyr::rename("chrom"=`#chrom`, "start"=chromStart,"end"=chromEnd) %>%
  GRanges()

# find overlaps between censat regions and methylation calls

chm13.ovls <- FindOvls(chm13_meth,censat.widths) %>%
  mutate(cell_line="chm13")
hg002.ovls <- FindOvls(hg002_meth,censat.widths) %>%
  mutate(cell_line="hg002") %>%
  GRanges()

# remove bad hg002 mapping regions
hg002.ovls <- as.data.frame(hg002.ovls[!hg002.ovls %over% hg002.flag,])

# violing plot of methylation frequency

censat.meth <- rbind(chm13.ovls,hg002.ovls)

ggplot(censat.meth, aes(x = methylated_frequency, y = name, fill=cell_line)) + geom_density_ridges(alpha=.5)+facet_wrap(~cell_line)+scale_fill_manual(values=c("hg002"="dodgerblue", "chm13"="magenta"))

ggsave(
  paste0(figs, "/", "AllSatv3_methylation_freq.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 10,
)
