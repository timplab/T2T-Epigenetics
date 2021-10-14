
## Script to make deeptools like heatmaps from methylation data in TEs

# loads libs 

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(GenomicRanges)
library(rtracklayer)
library(zoo)
# load meth 

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
chm13_meth <- import(paste0(dat, "/cut_and_run/210204_alignments/022021_CUTRUN/H3K27me3/CHM13_H3K27me3_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig"))

# what type of repeat 
type="SINE_HOTnot_"

# load rep file
reps <- read_tsv(paste0(dat, "/TE/HOT_TE/Deeptools_beds/", type, "CHM13ABpool_htmpSorted.bed")) %>%
  dplyr::rename("chr"=`#chrom`) %>% 
  filter(grepl("AluY",name)) %>%
  dplyr::rename("direction"=strand) %>%
  dplyr::select(c(chr, start, end, direction, deepTools_group)) %>%
  mutate(ID=row_number())

# set flanks 
flankn <- 3000
bodylen <- 1

# overlap with methylation calls and normalize by size 
l1_regions <- reps %>%
  mutate(start = start - flankn, end = end + flankn) %>%
  GRanges()

ovl <- findOverlaps(GRanges(chm13_meth), l1_regions)
genes.ovl <- as.data.frame(reps)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(GRanges(chm13_meth)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(direction == "+",start - gene_start, gene_end - start),
                dist = ifelse(dist < 0, dist/flankn,
                              ifelse(dist < genewidth,
                                     bodylen * dist / genewidth,
                                     bodylen + (dist - genewidth)/flankn)), 
                dist = round(dist,1)
  )


n_windows=500
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,deepTools_group) %>%
  summarise(med = median(score), top = quantile(score, 0.75), bot = quantile(score, 0.25), n_genes = length(score)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  group_by(deepTools_group) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 2, NA))


# plot frequency with loess --- if doing SINEs don't use loess unless u want to die 
freqplot <- ggplot(chm13.ovl.labs,aes(x=min,y=med_smooth,color=deepTools_group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))+geom_smooth(method="loess", span=.2)

# draw line where HOT elements end

# plot heatmap
plot <- ggplot(chm13.ovl,aes(x=dist,y=-ID,fill=score))+ geom_tile()+scale_fill_gradient(low = "blue", high = "red", na.value = NA, limits=c(0,1))+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))+theme(legend.position = "left", legend.direction="vertical")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())

# stack density and heatmap
plot_grid(freqplot, plot, ncol=1, align = "v", rel_heights=c(1/3,1))

ggsave(
  paste0(dat, "/figures/evol_meth/TE/",type, "_Heatmap.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 12,
)