#!/usr/bin Rscript

##### Script takes nanopore GpC methylation tabix indexed bed file, extracts and saves per read GpC methylation calls and bins GC methylation calls into 15 Kb bins ######################
library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(BSgenome.HG002.chrX)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/pooled/"

readsGC <- tabix_mbed(paste0(dat, "/HG002_nanonome_chrX_GpCmethylation_pooled.bed.gz"),extcol = "motif",by = "read")

size_sel <- readsGC %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 20000) %>%
  mutate(motif = "GC")

gccalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

write.csv(gccalls, paste0(dat, "/pooled/HG002_nanonome_chrX_GpCmethylation_ByCall.csv"))

calls.gr <-  gccalls %>%
  ungroup() %>%
  group_by(start) %>%
  summarise(num_meth = sum(mcall == 1), num_unmeth = sum(mcall == 0)) %>%
  mutate(chr = "chrX", start = start, end = start) %>%
  GRanges()


chrx.gr <- GRanges(seqinfo(BSgenome.HG002.chrX))


blocks <- genomeBlocks(BSgenome.HG002.chrX, chrs = seqnames(BSgenome.HG002.chrX), width = 15000)
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

score1 <- coverage(calls.gr, weight="num_meth")
score2 <- coverage(calls.gr, weight="num_unmeth")

binned_meth <- binnedSum(blocks, numvar = score1, "num_meth") %>%
  as.data.frame()
binned_unmeth <-binnedSum(blocks, numvar = score2, "num_unmeth")%>%
  as.data.frame()

meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  group_by(start, end) %>%
  mutate(sites = num_meth+num_unmeth) %>%
  mutate(freq = num_meth/sites) %>%
  na.omit()

freqmean=mean(meth_bins$freq)
freqsd=sd(meth_bins$freq)

gc_bins <- meth_bins %>%
  mutate(z_score = (freq - freqmean) / freqsd) %>%
  mutate(color = case_when(z_score > 0 ~ "Accessible", 
                           TRUE ~ "Inaccessible")) 

saveRDS(gc_bins, paste0(dat, "/pooled/ChrX_accessibilityZscore.rds"))