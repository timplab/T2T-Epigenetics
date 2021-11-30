#!/usr/bin/env Rscript 

# load libraries
library(tidyverse)
library(Repitools)
options(scipen=999)
library(zoo)
library(GenomicRanges)
library(BSgenome.t2t.v1.0.release)
library(optparse)

option_list = list(
  make_option(c("-i", "--nanopolish"), type="character", default=NULL, 
              help="path to nanopolish frequency  file"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="path to output directory"), 
  make_option(c("-w", "--window"), type="integer", default=NULL, 
              help="Smoothing window size"),
  make_option(c("-b", "--binsize"), type="integer", default=NULL, 
              help="Binning Size")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$nanopolish)){
  print_help(opt_parser)
  stop("provide path to nanopolish output.\n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Specify an output directory.\n", call.=FALSE)
}
if (is.null(opt$binsize)){
  print_help(opt_parser)
  stop("Specify a bin size.\n", call.=FALSE)
}

# input and output dir
dat=opt$out
win=opt$window
bin=opt$bin

# function for binned sum
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

# get chromosomes 
list=seqnames(BSgenome.t2t.v1.0.release)

# don't keep chrM
list=list[1:23]

# load nanopolish  output
chm13_meth <- read_tsv(opt$nanopolish) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome %in% list) %>%
  GRanges()

# split t2t genome into 10kb bins
blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = list, width = bin)


# find overlaps between all bed and 10kb bins, get bin average for all features 
chr.meth <- chm13_meth[seqnames(chm13_meth) %in% list]

score1 <- coverage(chr.meth, weight="called_sites_methylated")
score2 <- coverage(chr.meth, weight="called_sites_unmethylated")
score3 <- coverage(chr.meth, weight="num_motifs_in_group")


binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()

binned_unmeth <-binnedSum(blocks, numvar = score2, "called_sites_unmethylated")%>%
  as.data.frame()

binned_cov <- binnedSum(blocks, numvar = score3, "num_motifs_in_group") %>%
  as.data.frame()

# make meth bins and smooth with a rolling mean to make plot prettier
meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_cov, by = c("start", "end", "seqnames", "width", "strand")) %>%
  filter(seqnames %in% list) %>%
  # filter(start > rstart) %>%
  # filter(end < rend) %>%
  group_by(start, end, seqnames) %>%
  mutate(sites = called_sites_methylated+called_sites_unmethylated) %>%
  mutate(freq = called_sites_methylated/sites) %>%
  ungroup() %>%
  group_by(seqnames) %>%q
  arrange(start,seqnames) %>%
  mutate(smooth = rollmean(freq, win, fill = NA)) %>%
  ungroup() %>%
  arrange(seqnames, start)

write.table(meth_bins, paste0(opt$out, "/BinnedMethylationFrequency.tsv"), col.names = T, row.names = F, quote=F, sep = "\t")
