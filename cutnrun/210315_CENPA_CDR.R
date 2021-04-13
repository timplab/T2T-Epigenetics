
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
list=seqnames(BSgenome.t2t.v1.0.release)
# don't keep chrM, problems arise later because not all beds have chrM
list=list[1:23]

k=51 # 21 or 51
rep="022021" #022021 or 122920

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

name=paste0("_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_", k,"mers_single_mrg_meryl_sort.bigwig")
path=paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/210204_alignments/", rep, "_CUTRUN")

k51 <- read_tsv(paste0(dat, "/chm13_final_beds/chm13.v1.k", k,".single.mrg.markerdesert.bed"),col_names = c("chr", "start", "end", "name")) %>%
  mutate(num=1) %>%
  mutate(len=end-start) %>%
  mutate(end=start) %>%
  filter(chr %in% list) %>%
  GRanges()

hor = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621.bed"), col_names = F) %>%
  dplyr::filter(grepl("hor", X4))%>%
  mutate(status = ifelse(grepl("L", X4), "live", "dead")) %>%
  filter(status == "live") %>%
  group_by(X1) %>%
  rename(X2="start", X3="end") %>%
  rename(X1 = "chr") %>%
  GRanges() 


bw <- paste0(path,"/CENPA/CHM13_CA",name)
cenpa <- import(bw, format="bigwig")

bw <- paste0(path,"/IgG/CHM13_IgG",name)
igg <- import(bw, format="bigwig")

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated=called_sites-called_sites_methylated) %>%
  filter(chromosome %in% list) %>%
  GRanges()

hor.meth <- FindOvls(chm13_meth, hor) %>%
  GRanges()

score1 <- coverage(hor.meth, weight="called_sites_methylated")
score2 <- coverage(hor.meth, weight="called_sites_unmethylated")
score3 <- coverage(hor.meth, weight="num_motifs_in_group")
blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = list, width = 10000)

binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()
binned_unmeth <-binnedSum(blocks, numvar = score2, "called_sites_unmethylated")%>%
  as.data.frame()
binned_cov <- binnedMean(blocks, numvar = score3, "num_motifs_in_group") %>%
  as.data.frame()

meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_cov, by = c("start", "end", "seqnames", "width", "strand")) %>%
  filter(seqnames %in% list) %>%
  group_by(start, end, seqnames) %>%
  mutate(sites = called_sites_methylated+called_sites_unmethylated) %>%
  mutate(freq = called_sites_methylated/sites) %>%
  ungroup() %>%
  group_by(seqnames) %>%
  arrange(start,seqnames) %>%
  mutate(smooth = rollmean(freq, 3, fill = NA), site_smooth = rollmean(num_motifs_in_group, 3, fill = NA)) %>%
  ungroup() %>%
  arrange(seqnames, start)

# get the coordinates of the most hypomethylated region of every live HOR array
dips <- meth_bins %>%
  #filter(num_motifs_in_group > .01) %>%
  group_by(seqnames) %>%
  na.omit() %>%
  filter(freq == min(freq)) %>%
  mutate(dip=(end+start)/2) %>%
  select(c(seqnames,dip))

cenpa.hor <- FindOvls(cenpa,hor)

cenpa.hor.dist <- merge(cenpa.hor, dips, by="seqnames") %>%
  mutate(dist_cenpa = start-dip)


map.hor <- FindOvls(k51,hor)

map.hor.dist <- merge(map.hor, dips, by="seqnames") %>%
  mutate(dist_map = start-dip)


all.dist <- FindOvls(GRanges(cenpa.hor.dist),GRanges(map.hor.dist)) %>%
  mutate(dist=plyr::round_any(dist_cenpa, 10000)) %>%
  group_by(dist) %>%
  summarise(score=mean(score), len=mean(len))

breaks=c(1000,5000,10000,15000)
pal <- wes_palette("Zissou1", 5, type = "continuous")
ggplot()+geom_histogram(data=all.dist, aes(x=dist, y=score, fill=len),alpha=.8, stat="identity",size=2.5, position = "dodge" )+scale_fill_gradientn(colours = pal,name = paste0(k,"mer Desert Length"),breaks=breaks, limits=(c(1000,15000)))+theme(text = element_text(size=12))+labs(x="Distance from CDR", y = "CENPA Read Depth")+scale_x_continuous(breaks= c(-2000000,0,2000000), labels = c("-2 Mb","CDR","+2 Mb"))+xlim(-3000000,3000000)

ggsave(
  paste0(figs, "/", "CENPA_",k,"mer_metaplot_rep",rep, ".pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)

all.dist <- map.hor.dist %>%
  mutate(dist=plyr::round_any(dist_map, 20000)) %>%
  group_by(dist) %>%
  summarise(len=mean(len))


breaks=c(1000,5000,10000,15000)
pal <- wes_palette("Zissou1", 5, type = "continuous")
ggplot()+geom_histogram(data=all.dist, aes(x=dist, y=len),alpha=.8, stat="identity",size=2.5, position = "dodge")+labs(x="Distance from CDR", y = "Desert Length")+scale_x_continuous(breaks= c(-2000000,0,2000000), labels = c("-2 Mb","CDR","+2 Mb"))+xlim(-3000000,3000000)


ggsave(
  paste0(figs, "/", "CDR_desertLength",k,"mer_metaplot_rep",rep, ".pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)

igg.hor <- FindOvls(igg,hor)
igg.hor.dist <- merge(igg.hor, dips, by="seqnames") %>%
  mutate(dist_igg = start-dip)

all.dist.cenpa <- cenpa.hor.dist %>%
  mutate(dist=plyr::round_any(dist_cenpa, 20000)) %>%
  group_by(dist,seqnames) %>%
  summarise(score=mean(score))

all.dist.igg <- igg.hor.dist %>%
  mutate(dist=plyr::round_any(dist_igg, 20000)) %>%
  group_by(dist,seqnames) %>%
  summarise(score=mean(score))

ggplot()+geom_histogram(data=all.dist.cenpa, aes(x=dist, y=score),alpha=.8, stat="identity",size=2.5, position = "dodge",fill="red")+geom_histogram(data=all.dist.igg, aes(x=dist, y=score),alpha=.8, stat="identity",size=2.5, position = "dodge",fill="black")+theme(text = element_text(size=12))+labs(x="Distance from CDR", y = "CUT&RUN Read Depth")+scale_x_continuous(breaks= c(-2000000,0,2000000), labels = c("-2 Mb","CDR","+2 Mb"))+xlim(-3000000,3000000)

ggsave(
  paste0(figs, "/", "CDR_CENPA_IgG_",k,"mer_metaplot_rep",rep, ".pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)