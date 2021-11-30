library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(ggbreak)
options(scipen = 100)

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


hg002.gc <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_GpCmethylationFrequency.tsv")

hg002.flag <- hg002.gc %>%
  filter(chromosome != "chrY") %>%
  filter(called_sites < 10 | called_sites > 100) %>%
  mutate(start=start-500, end=start+500) %>%
  filter(end>start) %>%
  filter(start > 1) %>%
  dplyr::select(chromosome,start,end)


write.table(hg002.flag,
            "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Blacklist.bed", 
            sep="\t", quote=F, col.names = F, row.names = F)