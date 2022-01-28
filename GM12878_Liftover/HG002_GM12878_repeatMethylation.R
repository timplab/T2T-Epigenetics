source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
indir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/GM12878_liftover"
group1.freq.gm <- readRDS(paste0(indir, '/gm12878.methhap1.chm13.gr.rds')) %>%
  as.data.frame()  %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>%
  GRanges()
group2.freq.gm <- readRDS(paste0(indir, '/gm12878.methhap2.chm13.gr.rds')) %>%
  as.data.frame()  %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>%
  GRanges()

all.read.group1 <- readRDS("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/xci/HIGH_CGI_Reads.rds")
all.read.group2 <- readRDS("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/xci/LOW_CGI_Reads.rds")

group2.cg <- mbedByCall(all.read.group2) %>%
  drop_na(mcall)
group2.freq <- group2.cg %>%
  mutate(end=start) %>%
  group_by(chrom,start,end) %>%
  filter(n() > 5) %>%
  summarise(methylated_frequency=mean(mcall)) %>%
  GRanges()

group1.cg <- mbedByCall(all.read.group1) %>%
  drop_na(mcall)
group1.freq <- group1.cg %>%
  mutate(end=start) %>%
  group_by(chrom,start,end) %>%
  filter(n() > 5) %>%
  summarise(methylated_frequency=mean(mcall)) %>%
  GRanges()


l1s <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/xci/SINE_distance_to_XCI_genes.bed"), col_names = c("chr", "start", "end", "Rep", "score", "Direction", "start2", "end2", "score2", "num", "len", "start3", "HOT", "gene_chr", "gene_start", "gene_end","gene_strand", "name", "escape", "distance")) %>%
  filter(chr=="chrX") %>%
  mutate(ID=row_number()) %>%
  mutate(group_dist=case_when(distance==0 ~ 'geneic', 
                              abs(distance) > 0 & abs(distance) <= 5000 ~ "5kb", 
                              abs(distance)>5000 ~ "intergenic"))  %>%
  filter(group_dist != "intergenic")


l1s %>%
  group_by(HOT,group_dist) %>%
  summarise(sum=n()) %>%
  ungroup() %>%
  group_by(group_dist) %>%
  mutate(total=sum(sum)) %>%
  ungroup() %>%
  mutate(frac=sum/total) %>%
  ggplot(aes(x=group_dist, y=frac, fill=HOT))+geom_bar(stat="identity")

flankn <- 100
bodylen <- 500

# overlap with methylation calls and normalize by size 
l1_regions <- l1s %>%
  mutate(reg_end = ifelse(Direction == "+", end + flankn, start-flankn)) %>%
  mutate(reg_start = ifelse(Direction=="+", end - bodylen, start+bodylen)) %>%
  mutate(start=ifelse(Direction=="+", reg_start,reg_end)) %>%
  mutate(end=ifelse(Direction=="+", reg_end,reg_start)) %>%
  mutate(len=end-start) %>%
  GRanges()


ovl <- findOverlaps(GRanges(group1.freq), l1_regions)
genes.ovl <- as.data.frame(l1s)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) 

chm13.ovl1 <- as.data.frame(GRanges(group1.freq)[queryHits(ovl),])  %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",reg_end - start,start - reg_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )

ovl <- findOverlaps(GRanges(group2.freq), l1_regions)
genes.ovl <- as.data.frame(l1s)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) 

chm13.ovl2 <- as.data.frame(GRanges(group2.freq)[queryHits(ovl),])  %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",reg_end - start,start - reg_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )

chm13.ovl1 <- chm13.ovl1 %>%
  mutate(group="group1")
chm13.ovl2 <- chm13.ovl2  %>%
  mutate(group="group2")

chm13.ovl <- rbind(chm13.ovl1,chm13.ovl2)

n_windows=500
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group,escape,group_dist) %>%
  summarise(med = median(methylated_frequency),ID=ID) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  group_by(group,escape,group_dist) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 20, NA))

freqplot <- ggplot(chm13.ovl,aes(x=dist,y=methylated_frequency, color=group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+geom_smooth(method="loess",se=F, span=.3)+ylim(0,1)+facet_wrap(~escape+group_dist)

# plot heatmap
plot <- ggplot(chm13.ovl,aes(x=dist,y=-ID,fill=methylated_frequency))+ geom_tile()+scale_fill_gradient(low = "blue", high = "red", na.value = NA, limits=c(0,1))+theme(legend.position = "left", legend.direction="vertical")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())+facet_wrap(~escape)#+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))

# stack density and heatmap
plot_grid(freqplot, plot, ncol=1, align = "v", rel_heights=c(1/3,1))


stats <- chm13.ovl %>%
  dplyr::select(Rep, escape, group_dist) %>%
  distinct() %>%
  group_by(escape, group_dist) %>%
  summarise(n=n())
stats



flankn <- 100
bodylen <- 500

# overlap with methylation calls and normalize by size 
l1_regions <- l1s %>%
  mutate(reg_end = ifelse(Direction == "+", end + flankn, start-flankn)) %>%
  mutate(reg_start = ifelse(Direction=="+", end - bodylen, start+bodylen)) %>%
  mutate(start=ifelse(Direction=="+", reg_start,reg_end)) %>%
  mutate(end=ifelse(Direction=="+", reg_end,reg_start)) %>%
  mutate(len=end-start) %>%
  GRanges()


ovl <- findOverlaps(GRanges(group1.freq.gm), l1_regions)
genes.ovl <- as.data.frame(l1s)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) 

chm13.ovl1 <- as.data.frame(GRanges(group1.freq.gm)[queryHits(ovl),])  %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",reg_end - start,start - reg_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )

ovl <- findOverlaps(GRanges(group2.freq.gm), l1_regions)
genes.ovl <- as.data.frame(l1s)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) 

chm13.ovl2 <- as.data.frame(GRanges(group2.freq.gm)[queryHits(ovl),])  %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",reg_end - start,start - reg_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )

chm13.ovl1 <- chm13.ovl1 %>%
  mutate(group="group1")
chm13.ovl2 <- chm13.ovl2  %>%
  mutate(group="group2")

chm13.ovl <- rbind(chm13.ovl1,chm13.ovl2)

n_windows=500
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group,escape,group_dist) %>%
  summarise(med = median(methylated_frequency),ID=ID) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  group_by(group,escape,group_dist) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 20, NA))

freqplot <- ggplot(chm13.ovl,aes(x=dist,y=methylated_frequency, color=group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+geom_smooth(method="loess",se=F, span=.3)+ylim(0,1)+facet_wrap(~escape+group_dist)
freqplot
stats <- chm13.ovl %>%
  dplyr::select(Rep, escape, group_dist) %>%
  distinct() %>%
  group_by(escape, group_dist) %>%
  summarise(n=n())
stats




#l1s <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/xci/SINE_distance_to_XCI_genes.bed"), col_names = c("chr", "start", "end", "Rep", "score", #"Direction", "start2", "end2", "score2", "num", "len", "start3", "HOT", "gene_chr", "gene_start", "gene_end","gene_strand", "name", "escape", "distance")) %>%
#  filter(chr=="chrX") %>%
#  mutate(ID=row_number()) %>%
#  mutate(group_dist=case_when(distance==0 ~ 'geneic', 
#                              abs(distance) > 0 & abs(distance) <= 5000 ~ "5kb", 
#                              abs(distance)>5000 ~ "intergenic"))  %>%
#  filter(group_dist != "intergenic") %>%
#  GRanges()

rep="LTR"
l1s <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/xci/chm13_v1_RepeatMaskerV2_sorted_chrX_distance_to_XCI.bed"), col_names = c("chr", "start", "end", "Rep", "score", "Direction", "start2", "end2", "num", "len","rep_group", "rep_name", "divergence", "gene_chr", "gene_start", "gene_end","gene_strand", "name", "escape", "distance")) %>%
  filter(chr=="chrX") %>%
  filter(grepl(rep,rep_group)) %>%
  mutate(ID=row_number()) %>%
  mutate(group_dist=case_when(distance==0 ~ 'geneic', 
                              abs(distance) > 0 & abs(distance) <= 5000 ~ "5kb", 
                              abs(distance)>5000 ~ "intergenic"))  %>%
  filter(group_dist != "intergenic") %>%
  GRanges()


ovls.group1 <- FindOvls(GRanges(group1.freq), l1s) %>%
  mutate(hap="xa")
ovls.group2 <- FindOvls(GRanges(group2.freq), l1s)%>%
  mutate(hap="xi")

ovls.all <- rbind(ovls.group1,ovls.group2) %>%
  group_by(ID, group_dist,escape,hap) %>%
  summarize(mean_meth=mean(methylated_frequency))

p1 <- ggplot(ovls.all, aes(x=group_dist, y=mean_meth, fill=hap))+facet_wrap(~escape)+geom_boxplot()+stat_compare_means(method="wilcox.test")+ylim(0,1)

ovls.group1 <- FindOvls(GRanges(group1.freq.gm), l1s) %>%
  mutate(hap="xi")
ovls.group2 <- FindOvls(GRanges(group2.freq.gm), l1s)%>%
  mutate(hap="xa")

ovls.all.gm <- rbind(ovls.group1,ovls.group2) %>%
  group_by(ID, group_dist,escape,hap) %>%
  summarize(mean_meth=mean(methylated_frequency))

p2 <- ggplot(ovls.all.gm, aes(x=group_dist, y=mean_meth, fill=hap))+facet_wrap(~escape)+geom_boxplot()+stat_compare_means(method="wilcox.test")+ylim(0,1)

plot_grid(p1, p2, ncol=1)

ggsave(
  paste0(figs, "/", rep, "_methylation_CHM113_GM12878.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
)
