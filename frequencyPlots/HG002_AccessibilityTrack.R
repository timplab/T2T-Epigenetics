library(tidyverse)
library(rtracklayer)
library(BSgenome.t2t.v1.0.release)
library(Repitools)
library(GenomicRanges)

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

gcalls.gr <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_GpCmethylationFrequency.tsv")%>%
  filter(chromosome != "chrY") %>%
  GRanges()

win=200
chrs <- seqlevels(BSgenome.t2t.v1.0.release)

blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = chrs, width = win)

score1 <- coverage(gcalls.gr, weight="called_sites_methylated")
score2 <- coverage(gcalls.gr, weight="called_sites")

binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()

binned_all <-binnedSum(blocks, numvar = score2, "called_sites")%>%
  as.data.frame()

meth_bins <- merge(binned_meth, binned_all, by = c("start", "end", "seqnames", "width", "strand")) %>%
  group_by(seqnames,start, end) %>%
  mutate(freq = called_sites_methylated+called_sites) %>%
  na.omit()

freqmean=mean(meth_bins$freq)
freqsd=sd(meth_bins$freq)


gc_bins <- meth_bins %>%
  mutate(z_score = (freq - freqmean) / freqsd) %>%
  mutate(score=as.numeric(z_score)) %>%
  dplyr::select(-c(called_sites_methylated,called_sites,freq,z_score)) %>%
  GRanges()

seqinfo(gc_bins) <- seqinfo(BSgenome.t2t.v1.0.release)

export.bw(gc_bins, paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_", win,"Binned_Accessibility.bw"))
