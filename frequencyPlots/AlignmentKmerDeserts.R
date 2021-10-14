source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(GenomicRanges)
library(Repitools)
library(rtracklayer)

binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

cont <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/mappability/mappability_chm13v1_200bp.bw")

desert <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDesertsWithLen_50kbBinned.bed", col_names = c('chr', "start", "end", "num")) 

chm13_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/bam/ont.primary_meth50kb.bedgraph", col_names = c("chr", "start", "end", "cov")) %>%
  GRanges()

bsseq <-  read_tsv("/uru/Data/old_atium/Data/Nanopore/projects/ONT_HG002_bisulfite/bsseq/bismark/004_0111_001_R1_val_1_bismark_bt2_pe_sorted_50kavg.bed", col_names = c("seqnames", "start", "end", "cov")) %>%
  filter(seqnames != "chrY") %>%
  filter(seqnames != "chrX") %>%
  mutate(start=start+1)
bsseq$cov[ bsseq$cov == "." ] <- 0

desert$num[ desert$num == "." ] <- 0

desert <- desert %>%
  mutate(start=start+1) %>%
  mutate(num=as.numeric(num)) %>%
 # mutate(num=ifelse(num >= 50000, 50000, num)) %>%
  mutate(num=((num/50000)*100)) %>%
  GRanges()

chrs <- seqnames(BSgenome.t2t.v1.0.release)[1:23]


blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = chrs, width = 50000)
score <- coverage(chm13_meth, weight="cov")
binned_cov <- GenomicRanges::binnedAverage(blocks, numvar = score, "cov") %>%
  as.data.frame()

blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqnames(BSgenome.t2t.v1.0.release), width = 50000)

seqlevels(cont) <- seqlevels(blocks)
score1 <- coverage(cont, weight="score")


cont_binned <- GenomicRanges::binnedAverage(blocks, numvar = score1, "score") %>%
  as.data.frame()


dat <- merge(cont_binned,as.data.frame(desert), by=c("seqnames", "start", "end", "width")) %>%
  merge(binned_cov, by=c("seqnames", "start", "end", "width"))

dat2 <- merge(bsseq,as.data.frame(desert), by=c("seqnames", "start", "end")) 

p1 <- ggplot(dat, aes(x=num,y=score))+labs(y="Mappability Score", x="Percent Marker Desert")+
  xlim(0,100)+geom_point(alpha=.4, color="red", size=3)+
  theme(text = element_text(size=18))

p2 <- ggplot(dat, aes(x=num,y=cov))+labs(y="ONT Methylation Coverage", x="Percent Marker Desert")+
  xlim(0,100)+geom_point(alpha=.4, color="Blue", size=3)+
  theme(text = element_text(size=18))

p3 <- ggplot(dat2, aes(x=num,y=cov))+labs(y="Bisulfite Coverage", x="Percent Marker Desert")+
  xlim(0,100)+geom_point(alpha=.4, color="purple", size=3)+
  theme(text = element_text(size=18))+ylim(0,250)
plot_grid(p1, p2)

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures/DesertCompareBsseq.pdf",
  plot = p3,
  scale = 1,
  width = 5,
  height = 5,
)


ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures/DesertCompare.pdf",
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 5,
)

ont.deserts <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/ont.primary_meth50kb.MarkerDesertCoverage.bed", col_names = c("chr", "start","end","name","rchr","rstart","rend","cov")) %>%
  mutate(len=(end-start)) %>%
  group_by(name) %>%
  mutate(cov=max(cov)) %>%
  dplyr::select(-c(rchr,rstart,rend)) %>%
  distinct() %>%
  mutate(rlen=round(len,-3)) %>%
  group_by(rlen) %>%
  summarise(cov=median(cov)) %>%
  distinct() %>%
  ungroup() %>%
  group_by(rlen) %>%
  summarise(avgcov=mean(cov))

ggplot(ont.deserts, aes(x=rlen, y=avgcov))+geom_point()+scale_x_log10()+ylim(0,90)+geom_vline(xintercept=50000, linetype="dashed")+
  theme(text = element_text(size=18))+labs(x="Desert Length", y = "ONT Methylation Coverage")

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures/DesertONTcutoff.pdf",
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 5,
)

des <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDeserts.bed") %>%
  mutate(len=chromEnd-chromStart) %>%
  filter(len > 50000) 
