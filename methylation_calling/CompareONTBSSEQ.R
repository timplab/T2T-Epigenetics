library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg38)
t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)

chrs=seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:23]
hg38.cpg.loci <- findLoci(pattern = "CG",
                          subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                          include=chrs,
                          strand = "+")
length(hg38.cpg.loci)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"
hg002.nanopore <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv") %>%
  mutate(method="nanopore")

hg002.bsseq <- read_tsv("/uru/Data/old_atium/Data/Nanopore/projects/ONT_HG002_bisulfite/bsseq/bismark/CpGCalls/004_0111_001_R1_val_1_bismark_bt2_pe.bismark.cov.gz", col_names = c("chromosome", "start", "end", "methylated_frequency", "called_sites_methylated", "called_sites_unmethylated")) %>%
  mutate(start=start-1, end=end-1,methylated_frequency=methylated_frequency/100, called_sites =called_sites_methylated+called_sites_unmethylated, num_motifs_in_group=1, group_sequence="CG",method="bsseq") %>%
  dplyr::select(chromosome, start, end,num_motifs_in_group,called_sites,called_sites_methylated, methylated_frequency,group_sequence, method)

hg38_bsseq <- read_tsv("/uru/Data/old_atium/Data/Nanopore/projects/ONT_HG002_bisulfite/GRCh38_bismark/CpG.gz.bismark.cov.gz",col_names = c("chromosome", "start", "end", "methylated_frequency", "called_sites_methylated", "called_sites_unmethylated")) %>%
  mutate(start=start-1, end=end-1,methylated_frequency=methylated_frequency/100, called_sites =called_sites_methylated+called_sites_unmethylated, num_motifs_in_group=1, group_sequence="CG",method="bsseq_hg38") %>%
  dplyr::select(chromosome, start, end,num_motifs_in_group,called_sites,called_sites_methylated, methylated_frequency,group_sequence, method) %>%
  mutate(start=end) %>%
  filter(chromosome %in% chrs) %>%
  GRanges()

des <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDeserts.bed") %>%
  dplyr::rename('chr'=`#chrom`) %>%
  GRanges()
bsseq.gr <- GRanges(hg002.bsseq)
bsseq.ref <- FindOvls(t2t.cpg.loci,bsseq.gr) %>%
  select(-c(width,strand)) %>%
  dplyr::rename("chromosome"=seqnames)

grch38.ref <- FindOvls(hg38.cpg.loci,hg38_bsseq) %>%
  select(-c(width,strand)) %>%
  dplyr::rename("chromosome"=seqnames)

hg002.all <- rbind(hg002.nanopore,bsseq.ref,grch38.ref)

ggplot(hg002.all, aes(called_sites, fill = method, color=method, alpha=.5))+geom_histogram(position = 'identity', alpha=0.6, bins = 50)+labs(x="CpG Coverage", y="Count")+xlim(0,250)

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_bsseq_ONT_methCov.pdf",
  plot = last_plot(),
  width = 8,
  height = 5
)

hg002.low <- hg002.all %>%
  filter(called_sites <= 20)

ggplot(hg002.low, aes(called_sites, fill = method, color=method, alpha=.5))+geom_histogram(position = 'identity', alpha=0.6, bins = 20)+labs(x="CpG Coverage", y="Count")+xlim(0,20)

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_bsseq_ONT_methCovZoom.pdf",
  plot = last_plot(),
  width = 8,
  height = 5
)



total.stat <- hg002.all %>%
  group_by(method) %>%
  dplyr::summarize(total=sum(num_motifs_in_group)) %>%
  mutate(which="sampled")

method=c("bsseq","nanopore","bsseq_hg38")
total=as.numeric(c(length(t2t.cpg.loci),length(t2t.cpg.loci),length(hg38.cpg.loci)))
which=c("all","all","all")

all.cgs <- data.frame(method,total,which)
all.dat <- rbind(total.stat,all.cgs)
ggplot(all.dat, aes(x=method, y=total/1e6,fill=which, alpha=.5))+geom_bar(stat="identity", position="identity",fill=c("dodgerblue"))+theme_classic()

ggsave(
  paste0(figs,"/HG002_bsseq_ont_methbar.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8
)

hg002.nanopore.good <- hg002.nanopore %>%
  filter(called_sites > 20) %>%
  filter(called_sites <500)

hg002.bsseq.good <- hg002.bsseq %>%
  filter(called_sites > 20) %>%
  filter(called_sites < 500)

hg002.forplot <- merge(hg002.nanopore.good,hg002.bsseq.good, by = c("chromosome", "start", 'end')) 
  

hg002.forplot.gr <- GRanges(hg002.forplot)

des.meth <- FindOvls(hg002.forplot.gr,des)
nondes.meth <- hg002.forplot.gr[!(hg002.forplot.gr %over% des)]
nondes.meth.df <- as.data.frame(nondes.meth)

library(RColorBrewer)
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(nondes.meth.df$methylated_frequency.x, nondes.meth.df$methylated_frequency.y)
title <- sprintf("N = %d r = %.3f", nrow(nondes.meth.df), c)
ggplot(nondes.meth.df, aes(methylated_frequency.x, methylated_frequency.y)) +
  geom_bin2d(bins=30) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Bisulfite Methylation Frequency") +
  ylab("Nanopolish Methylation Frequency") +
  theme_bw(base_size=20) +
  ggtitle(title)

ggsave(filename = paste0(figs, "/HG002_bsseq_ONT_comparisonNonDesert.pdf"), plot = last_plot(),
       width = 5,
       height = 5)


library(RColorBrewer)
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(des.meth$methylated_frequency.x, des.meth$methylated_frequency.y)
title <- sprintf("N = %d r = %.3f", nrow(des.meth), c)
ggplot(des.meth, aes(methylated_frequency.x, methylated_frequency.y)) +
  geom_bin2d(bins=30) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Bisulfite Methylation Frequency") +
  ylab("Nanopolish Methylation Frequency") +
  theme_bw(base_size=20) +
  ggtitle(title)

ggsave(filename = paste0(figs, "/HG002_bsseq_ONT_comparisonDesert.pdf"), plot = last_plot(),
       width = 5,
       height = 5)
