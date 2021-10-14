library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(rtracklayer)
library(cowplot)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

h3k4me2 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/122920_CUTRUN/H3K4me2/CHM13_H3K4me2_cutnrun_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

h3k27me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/022021_CUTRUN/H3K27me3/CHM13_H3K27me3_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

#genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/NBPF.CHM13.combined.v4.MergedStranded.bed", col_names = c("chr", "start", "end", "strand", "name", "ID"))

genes <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF_duplicons.bed", col_names = c("chr", "start", "end"),delim=" ") %>%
  mutate(len=as.numeric(end)-as.numeric(start)) %>%
  mutate(ID=row_number()) %>%
  GRanges()

iso <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF.transcripts.meth.and.cutandrunMerged.bed", col_names = c("chr", "start", "end", "name", "score")) %>%
  dplyr::select(-c(chr,start,end))

iso <- merge(iso, genes, by="name")

tss <- genes %>%
  mutate(start = ifelse(strand == "+", start -500, end-500))%>%
  mutate(end = ifelse(strand == "+", start +500, end+500)) %>%
  GRanges()

meth.dat <- FindOvls(chm13_meth,genes) 
h3k4me2.dat <- FindOvls(h3k4me2,genes) 
h3k27me3.dat <- FindOvls(h3k27me3,genes)

meth.score <- meth.dat %>%
  group_by(seqnames,ID) %>%
  summarise(score=mean(methylated_frequency)) %>%
  mutate(mark="meth")%>%
  arrange(score)

cg.score <-  meth.dat %>%
  group_by(seqnames,ID) %>%
  summarise(score=n()/1000) %>%
  mutate(mark="CpG Density") %>%
  arrange(score)

h3k4me2.score <- h3k4me2.dat %>%
  group_by(seqnames,ID) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k4me2") 


h3k27me3.score <- h3k27me3.dat %>%
  group_by(seqnames,ID) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k27me3")

h3k27me3.score <- merge(h3k27me3.score,genes, by=c("ID"),all = TRUE) %>%
  dplyr::select(chr,start,end,name.y,ID,mark,score) %>%
  dplyr::mutate(score = replace_na(score, 0))

h3k4me2.score <- merge(h3k4me2.score,genes, by=c("ID"),all = TRUE) %>%
  dplyr::select(chr,start,end,name.y,ID,mark,score) %>%
  dplyr::mutate(score = replace_na(score, 0))

missing <- data.frame(ID=1:23,score=0)

ggplot(meth.score,(aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red")

p1 <- ggplot(meth.score,(aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red")

p2 <- ggplot(cg.score,(aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "black")

p3 <- ggplot(h3k4me2.score,(aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "orange")

p4 <- ggplot(h3k27me3.score,(aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple")

p5 <- ggplot(iso, (aes(x=1, y=as.factor(ID),fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "dodgerblue")

plot <- plot_grid(p1, p2,p3,p4, ncol=5)

ggsave(
  paste0(figs, "/NBPF_expressionForPhylo.pdf"),
  plot = plot,
  width = 20,
  height = 10
)
