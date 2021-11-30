library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(ggbreak)
options(scipen = 100)

hg002.nanopore <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled_50kb.tsv")  %>%
  GRanges()

#hg002.nanopore <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC.bw")

# load rep file
reps <- read_tsv(paste0(dat, "/revision_analysis/satellites/ACRO/annotate_units/ACRO_subunits_mergedStranded_withArray.bed"), col_names=c("chr", 'start', "end", "direction", "rep_chr", "rep_start", "rep_end", "type", "where","len", "array_status")) %>%
  filter(abs(end-start)>5500) %>%
  mutate(ID=row_number())

# set flanks 
flankn <- 50
bodylen <- 6200

# overlap with methylation calls and normalize by size 
l1_regions <- reps %>%
  filter(abs(end-start)>4000) %>%
  dplyr::mutate(gene_start = start, gene_end = end) %>%
  mutate(reg_end = ifelse(direction == "+", start+bodylen, end-bodylen)) %>%
  mutate(reg_start = ifelse(direction=="+", start-flankn, end+flankn)) %>%
  mutate(start=ifelse(direction=="+", reg_start,reg_end)) %>%
  mutate(end=ifelse(direction=="+", reg_end,reg_start)) %>%
  mutate(len=end-start) %>%
  GRanges()

ovl <- findOverlaps(GRanges(hg002.nanopore), l1_regions)
genes.ovl <- as.data.frame(reps)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(GRanges(hg002.nanopore)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist =  ifelse(direction == "+",start - gene_start, gene_end - start),
                dist = ifelse(dist < 0, dist/flankn,
                              ifelse(dist < genewidth,
                                     bodylen * dist / genewidth,
                                     bodylen + (dist - genewidth)/flankn)), 
                dist = round(dist,2)
  )

chm13.ovl %>%
  filter(array_status=="array") %>%
  ggplot(aes(x=dist,y=methylated_frequency, color=chr))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))+geom_smooth(method="loess",se=F, span=.1)+facet_wrap(~array_status)+ylim(0,1)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO"

ggsave(
  paste0(figs, "/", "ACRO_subunitMeta_CpG.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 5,
)
