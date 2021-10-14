library(tidyverse)
library(patchwork)
library(GenomicRanges)
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")

library(bsseq)

t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)

hg38.cpg.loci <- findLoci(pattern = "CG",
                          subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                          include=seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:23],
                            strand = "+")
length(hg38.cpg.loci)


hg38 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/chm13_to_grch38/methylation_frequency.tsv") %>%
  mutate(ref = "hg38")

chm13 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/all_readlengths/methylation_frequency.tsv")%>%
  mutate(ref = "chm13")

all.meth <- rbind(hg38, chm13)

all.meth <- all.meth %>%
  mutate(cov = called_sites/num_motifs_in_group)

plt1 <- ggplot(all.meth, aes(y=num_motifs_in_group,x=cov, fill = ref, color=ref, alpha=.5,stat="identity"))+geom_boxplot()+xlim(0,150)+labs(x="CpG Coverage", y="Count")+theme_classic()


plt2 <- ggplot(all.meth, aes(cov, fill = ref, color=ref, alpha=.5))+geom_histogram(position = 'identity', alpha=0.6, bins = 75)+xlim(0,150)+labs(x="CpG Coverage", y="Count")+theme_classic()

all.meth.low <- all.meth %>%
  mutate(cov = called_sites/num_motifs_in_group) %>%
  filter(cov <=20)


ggplot(all.meth.low, aes(cov, fill = ref, color=ref, alpha=.5))+geom_histogram(position = 'identity', alpha=0.6, bins = 20)+xlim(0,20)+labs(x="CpG Coverage", y="Count")+theme_classic()
ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/chm13_hg38_methCov0-20.pdf",
  plot = last_plot(),
  width = 8,
  height = 5
)

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/chm13_hg38_methCov.pdf",
  plot = plt2,
  width = 8,
  height = 5
)

total.stat <- all.meth %>%
  group_by(ref) %>%
  dplyr::summarize(total=sum(num_motifs_in_group)) %>%
  mutate(which="sampled")

ref=c("chm13","hg38")
total=as.numeric(c(length(t2t.cpg.loci),length(hg38.cpg.loci)))
which=c("all","all")

all.cgs <- data.frame(ref,total,which)

all.dat <- rbind(total.stat,all.cgs)
ggplot(all.dat, aes(x=ref, y=total/1e6,fill=which, alpha=.5))+geom_bar(stat="identity", position="identity",fill=c("dodgerblue"))+theme_classic()


ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/chm13_hg38_methbar.pdf",
  plot = last_plot(),
  width = 5,
  height = 8
)
