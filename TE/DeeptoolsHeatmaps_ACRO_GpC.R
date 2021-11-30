
## Script to make deeptools like heatmaps from methylation data in TEs

# loads libs 

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(GenomicRanges)
library(zoo)
# load meth 

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
#chm13_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
#  mutate(called_sites_unmethylated = called_sites - called_sites_methylated)

chm13_meth <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC.bw")

# load rep file
reps <- read_tsv(paste0(dat, "/revision_analysis/satellites/ACRO/annotate_units/ACRO_subunits_mergedStranded.bed"), col_names=c("chr", 'start', "end", "direction")) %>%
  dplyr::select(c(chr, start, end,direction)) %>%
  filter(abs(end-start)>4000) %>%
  mutate(ID=row_number())
# set flanks 
flankn <- 100
bodylen <- 6500

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

ovl <- findOverlaps(GRanges(chm13_meth), l1_regions)
genes.ovl <- as.data.frame(reps)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(GRanges(chm13_meth)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist =  ifelse(direction == "+",start - gene_start, gene_end - start),
                dist = ifelse(dist < 0, dist/flankn,
                              ifelse(dist < genewidth,
                                     bodylen * dist / genewidth,
                                     bodylen + (dist - genewidth)/flankn)), 
                dist = round(dist,2)
  )

n_windows=1000
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%

  group_by(cut) %>%
  summarise(med = median(score), top = quantile(score, 0.75), bot = quantile(score, 0.25), n_genes = length(score),ID=ID) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 30, NA))

# plot frequency with loess --- if doing SINEs don't use loess unless u want to die 
freqplot <- ggplot(chm13.ovl.labs,aes(x=min,y=med_smooth))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))+geom_smooth(method="loess",se=F, span=.1)
#+geom_line()#

# plot heatmap
plot <- ggplot(chm13.ovl,aes(x=dist,y=-ID,fill=score))+ geom_tile()+scale_fill_gradient(low = "white", high = "darkgreen", na.value = NA, limits=c(.2,.5))+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))+theme(legend.position = "left", legend.direction="vertical")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())+theme(legend.position = "bottom")

# stack density and heatmap
plot_grid(freqplot, plot, ncol=1, align = "v", rel_heights=c(1/3,1))

ggsave(
  paste0(dat, "/figures/evol_meth/TE/ACRO-HG002GpC_deeptoolsHeatmap.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 12,
)