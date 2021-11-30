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

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
figs=paste0(dat, "/figures")
t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)
rep <- read_tsv(paste0(dat, "/mike_annotation/chm13v1_repeat.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "Repeat")
sat <- read_tsv(paste0(dat, "/mike_annotation/chm13v1_satellite.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "Satellite")
sd <- read_tsv(paste0(dat, "/mike_annotation/chm13v1_segdup.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "SegDup")
line <- read_tsv(paste0(dat, "/mike_annotation/chm13v1_LINE.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "LINE/SINE")


rm <- rbind(rep, sat, sd, line) %>%
  GRanges()

chm13.ovls <- FindOvls(t2t.cpg.loci,rm) %>%
  group_by(rep_type) %>%
  summarize(chm13_CpGs=sum(width))

non_rep = length(t2t.cpg.loci) - sum(chm13.ovls$chm13_CpGs)
chm13.ovls <- chm13.ovls %>%
  add_row(rep_type="non-repetitive",chm13_CpGs=non_rep)

hg38.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                         include = seqlevels(BSgenome.t2t.v1.0.release),
                         strand = "+")
length(hg38.cpg.loci)

hg38.rep <- read_tsv(paste0(dat, "/mike_annotation/GRCh38p13_repeat (1).bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "Repeat")
hg38.sat <- read_tsv(paste0(dat, "/mike_annotation/GRCh38p13_satellite (1).bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "Satellite")
hg38.sd <- read_tsv(paste0(dat, "/mike_annotation/GRCh38p13_segdup.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "SegDup")
hg38.line <- read_tsv(paste0(dat, "/mike_annotation/GRCh38p13_LINE.bed"), col_names=F) %>%
  rename(X1="chr", X2="start", X3="end") %>%
  mutate(rep_type = "LINE/SINE")


hg38.rm <- rbind(hg38.rep, hg38.sat, hg38.sd, hg38.line) %>%
  GRanges()

hg38.ovls <- FindOvls(hg38.cpg.loci,hg38.rm) %>%
  group_by(rep_type) %>%
  summarize(hg38_CpGs=sum(width)) 

non_rep = length(hg38.cpg.loci) - sum(hg38.ovls$hg38_CpGs)

hg38.ovls <- hg38.ovls %>%
  add_row(rep_type="non-repetitive",hg38_CpGs=non_rep)

cpgs <- merge(chm13.ovls, hg38.ovls, by = "rep_type") %>%
  gather(asm, cpgs, chm13_CpGs:hg38_CpGs,-rep_type) %>%
  mutate(rep_type=if_else(cpgs < 60000, "other", rep_type)) %>%
  group_by(rep_type,asm) %>%
  summarise(cpgs=sum(cpgs))
  

repeatColors =c("LINE/SINE"="#7FC97F",
                "SegDup"="#BEAED4",
                "Satellite"="#FDC086", 
                "non-repetitive"="#FFFF99", 
                "Repeat"="#386CB0")

ggplot(cpgs, aes(x=asm, y=cpgs/1e6, fill= fct_reorder(rep_type,cpgs)))+geom_bar(stat="identity", color="black")+labs(x="Assembly", y="Number CpG sites (Millions)")+theme(text = element_text(size=20))+theme_classic()+scale_fill_manual("rep_type", values=repeatColors, drop=F)

ggsave(
  paste0(figs,"/CpG_numberPerRepeat_barV3.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 7
)

stat <- cpgs %>%
  group_by(asm) %>%
  mutate(total=sum(cpgs)) %>%
  group_by(asm, rep_type) %>%
  mutate(percent=(cpgs/total)*100)

bp<- ggplot(stat, aes(x="", y=percent, fill=rep_type))+
  geom_bar(width = 1, stat = "identity",color="black") +facet_wrap(~asm)
pie <- bp + coord_polar("y", start=0)+theme_void()+scale_fill_manual("rep_type", values=repeatColors, drop=F)
pie


stat.sum <- stat %>%
  select(-c(total, percent)) %>%
  spread(asm, cpgs) %>%
  mutate(dif=chm13_CpGs-hg38_CpGs) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(total=sum(dif)) 

bp<- ggplot(stat.sum, aes(x=rep_type, y=dif, fill=rep_type))+
  geom_bar(width = 1, stat = "identity",color="black")+theme_classic()+scale_fill_manual("rep_type", values=repeatColors, drop=F)

ggsave(
  paste0(figs,"/CpG_Difference_barV3.pdf"),
  plot = bp,
  scale = 1,
  width = 8,
  height = 5
)
