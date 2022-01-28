source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(rtracklayer)
library(Repitools)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

cenpa <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/022021_CUTRUN/CENPA/CHM13_CA_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrM") %>%
  GRanges()


# split t2t genome into 10kb bins
blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqnames(BSgenome.t2t.v1.0.release)[1:23], width = 10000)

score1 <- coverage(chm13_meth, weight="called_sites_methylated")
score2 <- coverage(chm13_meth, weight="called_sites_unmethylated")
score3 <- coverage(chm13_meth, weight="num_motifs_in_group")



binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()

binned_unmeth <-binnedSum(blocks, numvar = score2, "called_sites_unmethylated")%>%
  as.data.frame()

binned_cov <- binnedMean(blocks, numvar = score3, "num_motifs_in_group") %>%
  as.data.frame()

# make meth bins and smooth with a rolling mean to make plot prettier
meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_cov, by = c("start", "end", "seqnames", "width", "strand")) %>%
  # filter(start > rstart) %>%
  # filter(end < rend) %>%
  group_by(start, end, seqnames) %>%
  mutate(sites = called_sites_methylated+called_sites_unmethylated) %>%
  mutate(freq = called_sites_methylated/sites) %>%
  ungroup() %>%
  group_by(seqnames) %>%
  arrange(start,seqnames) %>%
  mutate(smooth = rollmean(freq, 3, fill = NA), site_smooth = rollmean(num_motifs_in_group, 3, fill = NA)) %>%
  ungroup() %>%
  GRanges()

cen <- read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.HOR_annotations.bed"), col_names = c("chr", "start", "end","name")) %>%
  mutate(status = ifelse(grepl("L", name), "live", "dead")) %>%
  filter(status=="live") %>%
  group_by(chr) %>%
  summarise(start=min(start), end=max(end)) %>%
  mutate(reg_end=end-start, reg_start=1) 

blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqnames(BSgenome.t2t.v1.0.release)[1:23], width = 10000)
seqlevels(cenpa) <- seqlevels(blocks)
score4 <- coverage(cenpa, weight="score")
binned_cenpa <- binnedSum(blocks, numvar = score4, "score") 

cenpa.cen <- FindOvls(binned_cenpa,GRanges(cen))%>%
  group_by(seqnames) %>%
  mutate(reg_start=min(start), reg_end=max(end)) %>%
  mutate(start=start-reg_start, end=start+width) 

meth.cen <- FindOvls(meth_bins,GRanges(cen))%>%
  group_by(seqnames) %>%
  mutate(reg_start=min(start), reg_end=max(end)) %>%
  mutate(start=start-reg_start, end=start+width) 

meth <- ggplot(meth.cen, aes(x = start/1e6, y=smooth))+geom_line() + labs(y="Methylation")+theme_classic(base_size = 10)+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)+facet_wrap(~seqnames,ncol = 1) + scale_fill_gradient(low="blue",high="red", limits=c(0,1))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+theme(legend.position = "bottom")


ggsave(
  paste0(figs, "/methyl_profiles/allChrMeth_PanelPlotv2.pdf"),
  plot = meth,
  scale = 1,
  width = 10,
  height = 15,
  useDingbats=FALSE
)

cenpa.plt <- ggplot(cenpa.cen, aes(xmin = start/1e6, xmax=end/1e6,ymin=0,ymax=1, fill=score))+geom_rect(size =1)+ labs(y="CENPA")+theme_classic(base_size = 10)+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+facet_wrap(~seqnames,ncol = 1)+ 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_fill_gradientn(colours = c("gray", "red", "darkred","black"), space = "rgb")

ggsave(
  paste0(figs, "/methyl_profiles/allChrCENPA_PanelPlotv2.pdf"),
  plot = cenpa.plt,
  scale = 1,
  width = 10,
  height = 8, 
  useDingbats=FALSE
)

