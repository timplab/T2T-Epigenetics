library(tidyverse)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/GM12878_liftover"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

hg002_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
  GRanges()
#chrX	76574877	76564829
gene_end=76574877+1000
gene_start=76564829

gene.freq <- as.data.frame(hg002_meth) %>%
  filter(seqnames=="chrX") %>%
  filter(start >= gene_start) %>%
  filter(end <= gene_end)


freqplot <- ggplot(gene.freq,aes(x=start,y=methylated_frequency))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+ylim(0,1)+geom_smooth(method="loess",se=F, span=.3)#+facet_wrap(~escape)#+geom_line()

ggsave(
  paste0(figs, "/HG002_TAF9B.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 5,
)




chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv"))


gene.freq <- as.data.frame(chm13_meth) %>%
  filter(chromosome=="chrX") %>%
  filter(start >= gene_start) %>%
  filter(end <= gene_end)

freqplot <- ggplot(gene.freq,aes(x=-start,y=methylated_frequency))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+ylim(0,1)+geom_smooth(method="loess",se=F, span=.3)#+facet_wrap(~escape)#+geom_line()

ggsave(
  paste0(figs, "/CHM13_PRKX.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 5,
)

