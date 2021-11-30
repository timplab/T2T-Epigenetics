#!/usr/bin/env Rscript

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()

hg002_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_CpGmethylationFrequency_20kb.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrY") %>%
  filter(chromosome != "chrM") %>%
  GRanges()

tar <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0_plus38Y_repeatmasker.out.bed"), col_names = c("chr", "start", "end", "rep")) %>%
  filter(rep == "TAR1") %>%
  filter(chr != "chrY") %>%
  GRanges()

cgi <- read_tsv("~/T2T-Epigenetics/CGI-t2t.txt") %>%
  filter(chr != "chrM") %>%
  GRanges()

ovls <- FindOvls(tar,cgi)

ends <- read_tsv(paste0(dat, "/reference/chm13.draft_v1.0.fasta.fai"), col_names = c("chr", "end")) %>%
  mutate(start = end - 20000) %>%
  mutate(side="start")

starts=data.frame(chr=ends$chr, start=rep(0, length(ends$chr))) %>%
  mutate(start=1, end=start + 20000) %>%
  mutate(side = "end")

telos <- rbind(ends, starts) %>%
  GRanges()

#write.table(as.data.frame(telos), paste0(dat, "/reference/telomere_coordinates.bed"), sep="\t", quote=F, row.names = F, col.names = F)
telocgis <- FindOvls(GRanges(ovls),telos) %>%
  mutate(ID=row_number())



chm13_meth_telo <- FindOvls(chm13_meth, GRanges(telocgis))%>%
  mutate(pos="telo")

chm13_avgmeth <- chm13_meth_telo %>%
  group_by(ID) %>%
  summarize(seqnames=seqnames,start=min(start), end=max(end), methylation=median(methylated_frequency),side=side) %>%
  distinct()
write.table(chm13_avgmeth, paste0(dat, "/CGI/telo_cgis/TeloCGI_meth.bed"), sep="\t", quote=F, row.names = F, col.names = T)
ggplot(chm13_meth_telo, aes(x=seqnames, y=methylated_frequency, fill=side))+geom_boxplot()

hg002_meth_telo <- FindOvls(hg002_meth, GRanges(telocgis))%>%
  mutate(pos="telo")

hg002_avgmeth <- hg002_meth_telo %>%
  group_by(ID) %>%
  summarize(seqnames=seqnames,start=min(start), end=max(end), methylation=median(methylated_frequency),side=side) %>%
  distinct()

write.table(hg002_avgmeth, paste0(dat, "/CGI/telo_cgis/HG002TeloCGI_meth.bed"), sep="\t", quote=F, row.names = F, col.names = T)

ggplot(hg002_meth_telo, aes(x=seqnames, y=methylated_frequency, fill=side))+geom_boxplot()

telocgis <- FindOvls(GRanges(ovls),telos)
#write.table(telocgis, paste0(dat, "/CGI/telo_cgis/TeloCGI.bed"), sep="\t", quote=F, row.names = F, col.names = T)
keepi <- findOverlaps(GRanges(ovls),telos)
freq.matched <- GRanges(ovls)[-queryHits(keepi)]

nontelo_cgi <- as.data.frame(freq.matched)
chm13_meth_nontelo <- FindOvls(chm13_meth, GRanges(nontelo_cgi))%>%
  mutate(pos="nontelo")%>%
  mutate(side=NA)
ggplot(chm13_meth_nontelo, aes(x=seqnames, y=methylated_frequency))+geom_boxplot()


hg002_meth_nontelo <- FindOvls(hg002_meth, GRanges(nontelo_cgi)) %>%
  mutate(pos="nontelo") %>%
  mutate(side=NA)
ggplot(hg002_meth_nontelo, aes(x=seqnames, y=methylated_frequency))+geom_boxplot()

hg002_all <- rbind(hg002_meth_telo, hg002_meth_nontelo)
ggplot(hg002_all, aes(x=pos, y=methylated_frequency))+geom_boxplot()

chm13_all <- rbind(chm13_meth_telo, chm13_meth_nontelo)
ggplot(chm13_all, aes(x=pos, y=methylated_frequency))+geom_boxplot()
