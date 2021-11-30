#!/usr/bin/env Rscript

### This script compares CpG site locations between GRCh38 and CHM13-T2T

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(cowplot)
library(tidyverse)
library(GenomicRanges)
library(bsseq)
library(BSgenome.t2t.v1.0.release)
library(BSgenome.Hsapiens.UCSC.hg38)

t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)

rm <- read_tsv(paste0(dat, "/annotations/chm13v1_polish-033121_minuschrM.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end", X4="rep_type", X5="len", X6="direction",X9="perc_div", X8="rep_class") %>%
  GRanges()

chm13.ovls <- FindOvls(t2t.cpg.loci,rm) %>%
  group_by(X7) %>%
  summarize(chm13_CpGs=sum(width))

non_rep = length(t2t.cpg.loci) - sum(chm13.ovls$chm13_CpGs)
chm13.ovls <- chm13.ovls %>%
  add_row(X7="non-repetitive",chm13_CpGs=non_rep)

hg38.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                         include = seqlevels(BSgenome.t2t.v1.0.release),
                         strand = "+")
length(hg38.cpg.loci)

hg38.rm <- read_tsv(paste0(dat, "/annotations/hg38_minuschrY-M_polish_033121.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end", X4="rep_type", X5="len", X6="direction",X9="perc_div", X8="rep_class") %>%
  GRanges()

hg38.ovls <- FindOvls(t2t.cpg.loci,hg38.rm) %>%
  group_by(X7) %>%
  summarize(hg38_CpGs=sum(width)) 

non_rep = length(hg38.cpg.loci) - sum(hg38.ovls$hg38_CpGs)

hg38.ovls <- hg38.ovls %>%
  add_row(X7="non-repetitive",hg38_CpGs=non_rep)

cpgs <- merge(chm13.ovls, hg38.ovls, by = "X7") %>%
  gather(asm, cpgs, chm13_CpGs:hg38_CpGs,-X7) %>%
  mutate(X7=if_else(cpgs < 60000, "other", X7)) %>%
  group_by(X7,asm) %>%
  summarise(cpgs=sum(cpgs))
  

repeatColors =c("DNA"="#C19936",
                "DNA?"="#C19935",
                "LINE"="#FFA500",
                "Low_complexity"="#75B043",
                "LTR"="#51B756",
                "RC"="#53BB74",
                "Retroposon"="#55BE9D",
                "RNA"="#ff4000",
                "rRNA"="#52BEBB",
                "scRNA"="#6B9AD3",
                "Simple_repeat"="#8992C9",
                "SINE"="pink",
                "snRNA"="#A886BC",
                "srpRNA"="#B67EB6",
                "Unspecified"="#C378B2",
                "tRNA"="#006400",
                "Unknown"="#BA55D3",
                "Satellite"="#53B0E4", 
                "non-repetitive"="dodgerblue")

ggplot(cpgs, aes(x=asm, y=cpgs/1e6, fill= fct_reorder(X7,cpgs)))+geom_bar(stat="identity", color="black")+scale_fill_manual("Repeat", values=repeatColors, drop=F)+labs(x="Assembly", y="Number CpG sites (Millions)")+theme(text = element_text(size=20))+theme_classic()

ggsave(
  paste0(figs,"/CpG_numberPerRepeat_barV2.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 7
)

stat <- cpgs %>%
  group_by(asm) %>%
  mutate(total=sum(cpgs)) %>%
  group_by(asm, X7) %>%
  mutate(percent=(cpgs/total)*100)

bp<- ggplot(stat, aes(x="", y=percent, fill=X7))+
  geom_bar(width = 1, stat = "identity",color="black") +facet_wrap(~asm)
pie <- bp + coord_polar("y", start=0)+scale_fill_manual("Repeat", values=repeatColors, drop=F)+theme_void()
pie


stat.sum <- stat %>%
  select(-c(total, percent)) %>%
  spread(asm, cpgs) %>%
  mutate(dif=chm13_CpGs-hg38_CpGs) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(total=sum(dif)) 

bp<- ggplot(stat.sum, aes(x=X7, y=dif, fill=X7))+
  geom_bar(width = 1, stat = "identity",color="black")+scale_fill_manual("Repeat", values=repeatColors, drop=F)+theme_classic()

ggsave(
  paste0(figs,"/CpG_Difference_barV2.pdf"),
  plot = bp,
  scale = 1,
  width = 8,
  height = 5
)
