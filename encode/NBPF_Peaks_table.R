library(tidyverse)

brain_h3k36me3 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_NBPF_intersections/NBPF_intersect_brain_microvascular_endothelial_cell_H3K36me3.bed", col_names=F) %>%
  mutate(sample="brain_h3k36me3")
brain_h3k27me3 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_NBPF_intersections/NBPF_intersect_brain_microvascular_endothelial_cell_H3K27me3.bed", col_names=F) %>%
  mutate(sample="brain_h3k27me3")

be2c_h3k27me3 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_NBPF_intersections/NBPF_intersect_BE2C_H3K27me3.bed", col_names=F) %>%
  mutate(sample="be2c_h3k27me3")
be2c_h3k36me3 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_NBPF_intersections/NBPF_intersect_BE2C_H3K36me3.bed", col_names=F) %>%
  mutate(sample="be2c_h3k36me3")

all.peaks <- rbind(brain_h3k36me3,brain_h3k27me3,be2c_h3k27me3,be2c_h3k36me3) %>%
  select(c(X4,sample)) %>%
  dplyr::rename("gene_name"=X4) %>%
  group_by(gene_name, sample) %>%
  summarise(n=n()) %>%
  spread(sample,n)  %>%
  ungroup %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) 
