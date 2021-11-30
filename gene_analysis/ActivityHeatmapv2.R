library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(rtracklayer)
library(cowplot)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()
hg002_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv") %>%
  GRanges()
hg002_gc <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC_zscore.bw") 


h3k4me2 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/122920_CUTRUN/H3K4me2/CHM13_H3K4me2_cutnrun_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

h3k27me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/022021_CUTRUN/H3K27me3/CHM13_H3K27me3_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF_Genes_Dups_intersect_isoseq.bed", col_names = c("chr", "gene_start", "gene_end", "strand", "gene_name", "dup_chr", "dup_start", "dup_end", "dup_strand", "dup_name", "ID", "iso")) 

tss <- genes %>%
  mutate(start=ifelse(strand=="+", gene_start-1000, gene_end-1000)) %>%
  mutate(end=ifelse(strand=="+", gene_start+1000, gene_end+1000)) %>%
  dplyr::select(c(chr,start,end,dup_name,iso,ID)) %>%
  GRanges()


hg002.gc.dat <- FindOvls(hg002_gc,tss) 
hg002.dat <- FindOvls(hg002_meth,tss) 
meth.dat <- FindOvls(chm13_meth,tss) 
h3k4me2.dat <- FindOvls(h3k4me2,tss) 
h3k27me3.dat <- FindOvls(h3k27me3,tss)

meth.score <- meth.dat %>%
  group_by(seqnames,ID,dup_name) %>%
  summarise(score=mean(methylated_frequency)) %>%
  mutate(mark="meth")%>%
  arrange(score)

hg002.gc.score <- hg002.gc.dat %>%
  group_by(seqnames,ID,dup_name) %>%
  summarise(score=mean(score)) %>%
  mutate(mark="hg002_gc")%>%
  arrange(score)

hg002.score <- hg002.dat %>%
  group_by(seqnames,ID,dup_name) %>%
  summarise(score=mean(methylated_frequency)) %>%
  mutate(mark="h002_meth")%>%
  arrange(score)


h3k4me2.score <- h3k4me2.dat %>%
  group_by(seqnames,ID,dup_name) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k4me2") 


h3k27me3.score <- h3k27me3.dat %>%
  group_by(seqnames,ID,dup_name) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k27me3")

p <- ggplot(hg002.gc.score,(aes(x=1, y=as.factor(ID*-1),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "pink")

p0 <- ggplot(hg002.score,(aes(x=1, y=as.factor(ID*-1),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red",limits=c(0,1))

p1 <- ggplot(meth.score,(aes(x=1, y=as.factor(ID*-1),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red",limits=c(0,1))

p3 <- ggplot(h3k4me2.score,(aes(x=1, y=as.factor(ID*-1),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "darkgreen")

p4 <- ggplot(h3k27me3.score,(aes(x=1, y=as.factor(ID*-1),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple")

p5 <- ggplot(as.data.frame(tss),(aes(x=1, y=as.factor(ID*-1),fill=iso)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "dodgerblue")

plot <- plot_grid(p, p0, p1,p3,p4,p5, ncol=5)
ggsave(
  paste0(figs, "/NBPF_expressionForPhyloTSSOnly.pdf"),
  plot = plot,
  width = 20,
  height = 10
)

