library(tidyverse)

peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/novel_genes/CHM13.novel.genes.v1.0.filtered.tss_activating_peaks.bed", col_names = c("chr", "start", "end", "type", "X5", "name", "start2", "end2", "num", "chr_2", "peak_start", "peak_end", "peak_name")) %>%
  separate(col=peak_name, into=c("cell_line", "mark"), sep="_H") %>%
  separate(col=mark, into=c("mark", "peak_ID"), sep="\\.") %>%
  mutate(mark = ifelse(grepl("3K27ac", mark), "H3K27ac", mark))%>%
  mutate(mark = ifelse(grepl("3K4me3", mark), "H3K4me3", mark)) %>%
  mutate(type=str_remove(type, "gene_biotype=")) %>%
  mutate(name=str_remove(name, "gene_name=")) %>%
  dplyr::select(c(name,type,mark,cell_line,peak_ID))

stats <- peaks %>%
  group_by(name,mark) %>%
  mutate(number_cell_line=length(unique(cell_line))) %>%
  filter(number_cell_line > 1)

stat.sum <- stats %>%
  select(c(name, number_cell_line, mark,type)) %>%
  distinct() %>%
  spread(mark, number_cell_line)%>% replace(is.na(.), 0)

type.stat <- stat.sum %>%
  ungroup() %>%
  select(c(type)) %>%
  group_by(type) %>%
  summarise(n=n())
  
