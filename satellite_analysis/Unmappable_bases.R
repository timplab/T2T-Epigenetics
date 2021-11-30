library(tidyverse)

nomap <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/mappability_chm13v1_200bp_score0.bed", col_names=c("chr", "start", "end", "score"))  %>%
  mutate(len=end-start)

total=sum(nomap$len)
