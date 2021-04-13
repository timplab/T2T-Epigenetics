
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
list=seqnames(BSgenome.t2t.v1.0.release)
# don't keep chrM, problems arise later because not all beds have chrM
list=list[1:23]


rep="122920" #022021 or 122920

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load cenpa data
k=51
name=paste0("_cutnrun_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_", k,"mers_single_mrg_meryl_sort.bigwig")
path=paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/210204_alignments/", rep, "_CUTRUN")

bw <- paste0(path,"/CENPA/CHM13_CA",name)
cenpa <- import(bw, format="bigwig")

bw <- paste0(path,"/IgG/CHM13_IgG",name)
igg <- import(bw, format="bigwig")


# get coordinates for just live HOR arrays
hor = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621.bed"), col_names = F) %>%
  dplyr::filter(grepl("hor", X4))%>%
  mutate(status = ifelse(grepl("L", X4), "live", "dead")) %>%
  filter(status == "live") %>%
  group_by(X1) %>%
  rename(X2="start", X3="end") %>%
  rename(X1 = "chr") %>%
  GRanges() 

cenpa.hor <- FindOvls(cenpa,hor)

# load kmer deserts
k=51
k51 <- read_tsv(paste0(dat, "/chm13_final_beds/chm13.v1.k", k,".single.mrg.markerdesert.bed"),col_names = c("chr", "start", "end", "name")) %>%
  mutate(num=1) %>%
  mutate(len=end-start) %>%
 # mutate(end=start) %>%
  filter(chr %in% list) %>%
  mutate(kmer="k51")

k=21
k21 <- read_tsv(paste0(dat, "/chm13_final_beds/chm13.v1.k", k,".single.mrg.markerdesert.bed"),col_names = c("chr", "start", "end", "name")) %>%
  mutate(num=1) %>%
  mutate(len=end-start) %>%
#  mutate(end=start) %>%
  filter(chr %in% list) %>%
  mutate(kmer="k21")

k100 <- read_tsv(paste0(dat, "/chm13_final_beds/chm13v1_unique_100mers_reformatted4overlaps.markerdesert.bed"),col_names = c("chr", "start", "end", "name")) %>%
  mutate(num=1) %>%
  mutate(len=end-start) %>%
#  mutate(end=start) %>%
  filter(chr %in% list) %>%
  mutate(kmer="k100")

# merge all deserts
kmers <- rbind(k51,k21,k100) %>%
  dplyr::rename("seqnames"=chr) %>%
  GRanges()

# bin cenpa data, 10kb bins with 100bp overlaps
score1 <- coverage(GRanges(cenpa.hor), weight="score")

sliding_blocks <- slidingWindows(hor, width = 10000, step = 100L)
sliding_blocks <- unlist(sliding_blocks)

# find biggest cenpa peak, set peak coord to middle of bin
binned_cenpa <- binnedSum(sliding_blocks, numvar = score1, "score") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(peak=(end+start)/2) %>%
  group_by(seqnames) %>%
  filter(score == max(score)) %>%
  mutate(peak =(start+end)/2) %>%
  summarise(peak=mean(peak))

# pull out kmers that are in live hors
kmers.hor <- FindOvls(kmers,hor)

# find distance of kmer desert to peak
cenpa.hor.dist <- merge(kmers.hor, binned_cenpa, by="seqnames") %>%
  group_by(seqnames) %>%
  mutate(dist = start-peak) %>%
  mutate(cordist=ifelse(dist > 0, start-peak, end-peak)) %>%
  group_by(kmer) %>%
  mutate(dist_round=plyr::round_any(cordist, 50000)) %>%
  group_by(dist_round,kmer,seqnames) %>%
  summarise(len=max(len))

# reorder levels for plot
cenpa.hor.dist$kmer <- factor(cenpa.hor.dist$kmer, levels=c("k21", "k51", "k100"))

pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot()+geom_histogram(data=cenpa.hor.dist, aes(x=dist_round, y=len, fill=kmer), alpha=1, stat="identity", position="dodge",size=2.5)+theme(text = element_text(size=12))+labs(x="Distance from CENPA Peak", y = "Length of kmer Desert")+scale_x_continuous(breaks= c(-1000000,0,1000000), labels = c("-1 Mb","CENPA","+1 Mb"))+xlim(-2000000,2000000)+facet_wrap(~kmer)


ggsave(
  paste0(figs, "/", "CENPA_allKmer_metaplot_rep",rep, ".pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 5,
)


# with read depth on y-axis just for k51
cenpa.hor.dist <- merge(cenpa.hor, binned_cenpa, by="seqnames") %>%
  mutate(dist_cenpa = start-peak)

map.hor <- FindOvls(GRanges(k51),hor)


all.dist <- FindOvls(GRanges(cenpa.hor.dist),GRanges(map.hor)) %>%
  group_by(seqnames) %>%
  mutate(dist=plyr::round_any(dist_cenpa, 20000)) %>%
  group_by(dist) %>%
  summarise(score=mean(score), len=max(len))

pal <- wes_palette("Zissou1", 50, type = "continuous")

ggplot()+geom_histogram(data=all.dist, aes(x=dist, y=score, fill=len),alpha=1, stat="identity",size=2.5)+theme(text = element_text(size=12))+labs(x="Distance from CENPA Peak", y = "CENPA Read Depth")+scale_x_continuous(breaks= c(-1000000,0,1000000), labels = c("-1 Mb","CENPA","+1 Mb"))+scale_fill_gradientn(colours = pal, name="Max 51mer Desert Size")+xlim(-2000000,2000000)

ggsave(
  paste0(figs, "/", "CENPAonly_51Kmer_metaplot_rep",rep, ".pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 5,
)
