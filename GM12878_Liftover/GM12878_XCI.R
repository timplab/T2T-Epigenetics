source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(zoo)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/GM12878_liftover"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


protein_coding <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/mitchell_merged_bed/protein_coding_names.txt"),col_names="gene") 


regs.all <- read_delim(paste0(dat, "/revision_analysis/promoter_clustering/mitchell_merged_bed/chm13_v1_tss.CGI_ntranscripts_intersections.bed"), col_names=c("chr","start","end","strand","name","gene_start","gene_end", "transcripts", "cgi_chr","cgi_start","cgi_end"),delim="\t")  %>%
  filter(name %in% protein_coding$gene) %>%
  mutate(quantile=case_when(transcripts == 0 ~ "no_expression", 
                            transcripts > 0 & transcripts < 113 ~ "medium",
                            transcripts >=100 ~ "high")) %>%
  filter(chr=="chrX") 

regs.cgi <- regs.all  %>%
  mutate(region="cgi") 

regs.gb <- regs.all  %>%
  mutate(start=ifelse(strand=="+", cgi_end,gene_end))%>%
  mutate(end=ifelse(strand=="+", gene_end,cgi_end)) %>%
  mutate(len=end-start) %>%
  filter(chr=="chrX") %>%
  filter(len > 100) %>%
  dplyr::select(-c(len))  %>%
  mutate(region="gb")


regs <- rbind(regs.gb,regs.cgi) %>%
  distinct() %>%
  # mutate(start=ifelse(strand=="+", gene_start,gene_end))%>%
  #  mutate(end=ifelse(strand=="+", gene_end,gene_start)) %>%
  GRanges()


known_escapers <- read_csv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/chrX/XCI_2014_paper.csv", col_names = c('name', 'escape')) %>%
  filter(escape %in% c("escape from XCI in 27 tissues","subject to XCI in 27 tissues"))


flankn=1000
bodylen <- 4
# add methylation data
l1_regions <- as.data.frame(regs) %>%
  mutate(start = start - flankn, end = end + flankn) %>%
  GRanges()


group1.freq <- readRDS(paste0(figs, '/gm12878.methhap1.chm13.gr.rds'))
group2.freq <- readRDS(paste0(figs, '/gm12878.methhap2.chm13.gr.rds'))


ovl <- findOverlaps(group1.freq, GRanges(l1_regions))
genes.ovl <- as.data.frame(regs)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) %>%
  mutate("Direction"=strand)

# calculate distance from start and end, depending on strand 
chm13.ovl.group1 <- as.data.frame(group1.freq[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(Direction=="+", start - reg_start, reg_end-start)) %>%
  mutate(dist = ifelse(dist < 0, dist/flankn,
                       ifelse(dist < genewidth,
                              bodylen * dist / genewidth,
                              bodylen + (dist - genewidth)/flankn)), 
         dist = round(dist,2)
  ) %>%
  mutate(group = "group1")




ovl <- findOverlaps(group2.freq, GRanges(l1_regions))
genes.ovl <- as.data.frame(regs)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(reg_start = start, reg_end = end) %>%
  mutate("Direction"=strand)

# calculate distance from start and end, depending on strand 
chm13.ovl.group2 <- as.data.frame(group2.freq[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(Direction=="+", start - reg_start, reg_end-start)) %>%
  mutate(dist = ifelse(dist < 0, dist/flankn,
                       ifelse(dist < genewidth,
                              bodylen * dist / genewidth,
                              bodylen + (dist - genewidth)/flankn)), 
         dist = round(dist,2)
  ) %>%
  mutate(group = "group2")

chm13.ovl <- rbind(chm13.ovl.group1,chm13.ovl.group2) %>%
  merge(known_escapers, by="name")

#write.table(chm13.ovl, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis#/promoter_clustering/ClusteredMethylationFrequency_chrX_proteinCoding_withUncertainty.tsv",
#            quote=F,
#            row.names=F,
#            sep="\t")
#
# bin data to make plot pretty
n_windows=5000
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group,escape) %>%
  summarise(med = median(methylated_frequency), top = quantile(methylated_frequency, 0.75), bot = quantile(methylated_frequency, 0.25), n_genes = length(methylated_frequency)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  group_by(group,escape) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 15, NA),top_smooth = rollmean(top, 10, NA),bot_smooth = rollmean(bot, 10, NA))

# plot -- include annotation for number of repeats in each group
p <- ggplot(chm13.ovl.labs,aes( x = min, y = med_smooth, color=group, fill=group), alpha=.5)+
  geom_line(aes(y=med_smooth), alpha=.5, size=1) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = bodylen)+
  labs( x = "Genomic Position", y = "Aggregated Methylation Frequency") +
  theme_classic()+facet_wrap(~escape,ncol = 1)+ylim(0,1)

#+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb")))+theme(legend.position = 'none')+facet_wrap(~escape,ncol = 1)+ylim(0,1)
p


ggsave(
  paste0(figs, "/GM12878_chrX_metaplot.pdf"),
  plot = p,
  scale = 1,
  width = 12,
  height = 8,
)




rep="L1HS"

l1s <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/revision/new_deeptools_beds/LINE_HOTnot_CHM13ABpool_anti-sense_htmpSorted_3Prime.bed") %>%
  dplyr::rename("chr"=`#chrom`) %>%
  dplyr::rename("Direction"=strand) %>%
  mutate(len=abs(end-start))
gene_regs <- as.data.frame(regs) %>%
  mutate(start=ifelse(strand=="+", gene_start,gene_end)) %>%
  mutate(end=ifelse(strand=="+", gene_end,gene_start)) %>%
  merge(known_escapers, by="name")


l1_genes <- FindOvls(GRanges(l1s), GRanges(gene_regs)) %>%
  dplyr::select(-c(gene_start, gene_end)) %>%
  mutate(ID=row_number())

flankn <- 100
bodylen <- 7000

# overlap with m
# overlap with methylation calls and normalize by size 
l1_regions <- l1_genes %>%
  dplyr::mutate(gene_start = start, gene_end = end) %>%
  mutate(reg_end = ifelse(Direction == "+", end + flankn, start-flankn)) %>%
  mutate(reg_start = ifelse(Direction=="+", end - bodylen, start+bodylen)) %>%
  mutate(start=ifelse(Direction=="+", reg_start,reg_end)) %>%
  mutate(end=ifelse(Direction=="+", reg_end,reg_start)) %>%
  GRanges()

ovl <- findOverlaps(GRanges(group1.freq), l1_regions)
genes.ovl <- as.data.frame(l1_genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl1 <- as.data.frame(GRanges(group1.freq)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",gene_end - start,start - gene_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )


ovl <- findOverlaps(GRanges(group2.freq), l1_regions)
genes.ovl <- as.data.frame(l1_genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl2 <- as.data.frame(GRanges(group2.freq)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",gene_end - start,start - gene_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,10)*-1)
#  )

chm13.ovl1 <- chm13.ovl1 %>%
  mutate(group="group1")
chm13.ovl2 <- chm13.ovl2  %>%
  mutate(group="group2")

chm13.ovl <- rbind(chm13.ovl1,chm13.ovl2) %>%
  group_by(group) %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated))

n_windows=1000
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group,escape) %>%
  summarise(med = median(methylated_frequency),ID=ID,escape=escape) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  group_by(group,escape) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 2, NA))

# plot frequency with loess --- if doing SINEs don't use loess unless u want to die 
freqplot <- ggplot(chm13.ovl.labs,aes(x=min,y=med_smooth, color=group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+ylim(0,1)+ geom_smooth(method="loess",se=F, span=.3)+facet_wrap(~escape)#+geom_line()

# plot heatmap
plot <- ggplot(chm13.ovl,aes(x=dist,y=-ID,fill=methylated_frequency))+ geom_tile()+scale_fill_gradient(low = "blue", high = "red", na.value = NA, limits=c(0,1))+theme(legend.position = "left", legend.direction="vertical")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())+facet_wrap(~group+escape,ncol=1)#+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))

# stack density and heatmap
plot_grid(freqplot, plot, ncol=1, align = "v", rel_heights=c(1/3,1))

ggsave(
  paste0(figs, "/", rep, "_Heatmap.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 12,
)

group1 <- group1.freq %>%
  as.data.frame() %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>% 
  mutate(group="hap1")
group2 <- group2.freq %>%
  as.data.frame() %>%
  mutate(methylated_frequency=Methylated/(Methylated+Unmethylated)) %>% 
  mutate(group="hap2")

#chrX:23284669-23326885
#3333049	3218669
#TAF9B chrX:76564830-76574732
methylation_freq <- rbind(group1,group2)
gene_end=3333049+1000
gene_start=3300000

gene.freq <- methylation_freq %>%
  filter(seqnames=="chrX") %>%
  filter(start >= gene_start) %>%
  filter(end <= gene_end)


n_windows=1000
gene.freq$cut = cut(gene.freq$start, breaks=n_windows)

chm13.ovl.labs <- gene.freq %>%
  group_by(cut,group) %>%
  summarise(med = median(methylated_frequency)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  group_by(group) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 2, NA))

freqplot <- ggplot(gene.freq,aes(x=start,y=methylated_frequency, color=group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+ylim(0,1)+geom_smooth(method="loess",se=F, span=.3)#+facet_wrap(~escape)#+geom_line()

ggsave(
  paste0(figs, "/GM12878_PRKX.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 5,
)
