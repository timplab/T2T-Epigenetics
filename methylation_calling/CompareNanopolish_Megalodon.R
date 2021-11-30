library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"
chm13.nanopolish<-  read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  filter(called_sites >10) %>%
  dplyr::select(chromosome,start,methylated_frequency) 

chm13.meg<-  read_tsv(paste0(dat, "/megalodon/megalodon_CpG/modified_bases.re_aggregated.5mC.bed"), col_names=c("chromosome", "start", "end", "score", "cov","strand", "start1", "end1", "nums", "cov", "methylated_frequency"))

chm13.meg.merge <- chm13.meg %>%
  mutate(start=ifelse(strand=="+", start, start-1)) %>%
  mutate(end=ifelse(strand=="+", end, end-1)) %>%
  group_by(chromosome,start,end) %>%
  summarise(methylated_frequency=mean(methylated_frequency), cov=sum(cov)) %>%
  filter(cov>10)%>%
  dplyr::select(chromosome,start,methylated_frequency) 

chm13.all <- merge(chm13.meg.merge,chm13.nanopolish, by=c("chromosome","start"))
library(RColorBrewer)
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(chm13.all$methylated_frequency.x, chm13.all$methylated_frequency.y)
title <- sprintf("N = %d r = %.3f", nrow(chm13.all), c)
ggplot(chm13.all, aes(methylated_frequency.x/100, methylated_frequency.y)) +
  geom_bin2d(bins=30) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Megalodon") +
  ylab("Nanopolish") +
  theme_bw(base_size=20) +
  ggtitle(title)

ggsave(filename = paste0(figs, "/CHM13_Nanopolish_Megalodon.pdf"), plot = last_plot(),
       width = 5,
       height = 5)