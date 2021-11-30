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
hg002_gc <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_GpCmethylationFrequency_20kb.bigwig") 


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

tss <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/whole_genome/NBPF_promoter_CGIs.bed", col_names=F) %>%
  select(X6, X7, X8, X9) %>%
  distinct() %>%
  dplyr::rename("gene"=X6, "chr"=X7, "start"=X8, "end"=X9) %>%
  mutate(ID=row_number()) 
nbpf8 <- data.frame(gene="NBPF8", chr="chr1", start=146655204,end=146656214, ID=10)

tss <- rbind(tss,nbpf8) %>%
  GRanges()

be2c_h3k36me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/bam_pileups/BE2C_H3K36me3.chm13v1.bw")
be2c_h3k27me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/bam_pileups/BE2C_H3K27me3.chm13v1.bw")

brain_microvascular_endothelial_cell_h3k36me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/bam_pileups/brain_microvascular_endothelial_cell_H3K36me3.chm13v1.bw")
brain_microvascular_endothelial_cell_h3k27me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/bam_pileups/brain_microvascular_endothelial_cell_H3K27me3.chm13v1.bw")

hg002.gc.dat <- FindOvls(hg002_gc,tss) 
hg002.dat <- FindOvls(hg002_meth,tss) 
meth.dat <- FindOvls(chm13_meth,tss) 
h3k4me2.dat <- FindOvls(h3k4me2,tss) 
h3k27me3.dat <- FindOvls(h3k27me3,tss)
be2c_h3k36me3.dat <- FindOvls(be2c_h3k36me3,tss)
be2c_h3k27me3.dat <- FindOvls(be2c_h3k27me3,tss)
brain_h3k36me3.dat <- FindOvls(brain_microvascular_endothelial_cell_h3k36me3,tss)
brain_h3k27me3.dat <- FindOvls(brain_microvascular_endothelial_cell_h3k27me3,tss)

be2c_h3k36me3.score <- be2c_h3k36me3.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="be2c_h3k36me3") 
be2c_h3k27me3.score <- be2c_h3k27me3.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="be2c_h3k27me3") 
brain_h3k27me3.score <- brain_h3k27me3.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="brain_h3k27me3") 
brain_h3k36me3.score <- brain_h3k36me3.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="brain_h3k36me3") 

meth.score <- meth.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=mean(methylated_frequency)) %>%
  mutate(mark="meth")%>%
  arrange(score)

hg002.gc.score <- hg002.gc.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=mean(score)) %>%
  mutate(mark="hg002_gc")%>%
  arrange(score)

hg002.score <- hg002.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=mean(methylated_frequency)) %>%
  mutate(mark="h002_meth")%>%
  arrange(score)

cg.score <-  meth.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=n()/1000) %>%
  mutate(mark="CpG Density") %>%
  arrange(score)

h3k4me2.score <- h3k4me2.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k4me2") 


h3k27me3.score <- h3k27me3.dat %>%
  group_by(seqnames,ID,gene) %>%
  summarise(score=max(score)) %>%
  mutate(mark="h3k27me3")

p <- ggplot(hg002.gc.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "pink", limits=c(0,.4))

p0 <- ggplot(hg002.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,1))

p1 <- ggplot(meth.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,1))

p2 <- ggplot(cg.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "black")

p3 <- ggplot(h3k4me2.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "darkgreen", limits=c(0,500))

p4 <- ggplot(h3k27me3.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple", limits=c(0,200))

p5 <- iso %>%
  filter(name %in% as.data.frame(tss)$gene)%>%
  ggplot((aes(x=1, y=name,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "dodgerblue")


p6 <- ggplot(be2c_h3k36me3.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "orange", limits=c(0,100))

p7 <- ggplot(be2c_h3k27me3.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple", limits=c(0,200))
p8 <- ggplot(brain_h3k36me3.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "orange", limits=c(0,100))

p9 <- ggplot(brain_h3k27me3.score,(aes(x=1, y=gene,fill=score)))+geom_tile()+
  scale_fill_gradient(low = "white", high = "purple", limits=c(0,200))

  
plot <- plot_grid(p, p0, p1,p3,p4,p5, ncol=5)
plot2 <- plot_grid(p6, p7,p8,p9, ncol=4)
ggsave(
  paste0(figs, "/NBPF_expressionForPhyloCGIOnly.pdf"),
  plot = plot,
  width = 20,
  height = 10
)

ggsave(
  paste0(figs, "/NBPF_expressionForPhyloCGIOnlyENCODE.pdf"),
  plot = plot2,
  width = 20,
  height = 10
)