library(tidyverse)
library(rtracklayer)
library(BSgenome.t2t.v1.0.release)

info = seqinfo(BSgenome.t2t.v1.0.release)[c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8", "chr9","chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22","chrX")]

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


dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


hg002 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_CpGmethylationFrequency.tsv")

hg002.flag <- hg002 %>%
  filter(chromosome != "chrY") %>%
  filter(called_sites < 10 | called_sites > 100) %>%
  GRanges()

chm13.flag <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  filter(chromosome != "chrM") %>%
  filter(called_sites < 10 | called_sites > 100) %>%
  GRanges()

chrs <- seqnames(BSgenome.t2t.v1.0.release)[1:23]
blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = chrs, width = 50000)

score <- coverage(hg002.flag, weight="num_motifs_in_group")
binned_hg002 <- binnedSum(blocks, numvar = score, "num_motifs_in_group") %>%
  as.data.frame() %>%
  na.omit() %>%
  dplyr::rename("score"=num_motifs_in_group) %>%
  mutate(score=as.numeric(score)) %>%
  GRanges()

score2 <- coverage(chm13.flag, weight="num_motifs_in_group")
binned_chm13 <- binnedSum(blocks, numvar = score2, "num_motifs_in_group") %>%
  as.data.frame() %>%
  na.omit() %>%
  dplyr::rename("score"=num_motifs_in_group) %>%
  mutate(score=as.numeric(score)) %>%
  GRanges()

seqinfo(binned_hg002) <- info
seqinfo(binned_chm13) <- info
export.bw(binned_hg002, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_50kbBinned_flaggedCalls.bw")

export.bw(binned_chm13, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/CHM13_50kbBinned_flaggedCalls.bw")
