library(tidyverse)
library(bsseq)
library(biovizBase)
library(Repitools)
library(Biostrings)
library(ggplot2)
library(BSgenome)
library(zoo)
library(GenomicRanges)
library(BSgenome.t2t.v1.0.release)

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

###functions#####

CalculateCpG <- function (obj, ..., step, width, as.prob, as.array, fast.moving.side,with.labels) 
{
  require(BSgenome)
  seqs <- getSeq(obj, ..., as.character = FALSE)
  if (missing(step)) {
    res <- oligonucleotideFrequency(seqs, 2, as.prob = as.prob)[,7]
  }
  else {
    res <- lapply(seqs, function(x) oligonucleotideFrequency(seqs, width = 2, step = step, as.prob = as.prob, as.array = as.array,fast.moving.side=fast.moving.side,with.labels=with.labels)[,7])
    if (length(res) == 1) 
      res <- unlist(res)
  }
  res
}
########


# build t2t BSgenome package and load it as BSgenome.t2t.v1.1

# split genome into 200bp bins and calculate CpG density per bin, save as GRanges object for easy reloading, run once comment out

chm13.gr <- GRanges(seqinfo(BSgenome.t2t.v1.0.release))
blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqnames(BSgenome.t2t.v1.0.release), width = 200)

sliding_blocks <- slidingWindows(blocks, width = 199, step = 1L)
sliding_blocks <- unlist(sliding_blocks)

chm13_CG <- CalculateCpG(BSgenome.t2t.v1.0.release, sliding_blocks, as.prob = F)

chm13 <- GRanges(sliding_blocks, CpG = chm13_CG)


#saveRDS(chm13, file =paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))
chm13 <- readRDS(paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))

freq = read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb.tsv"))
#bis = read.bismark(paste0(dat, "/methylation_calls/CpGmethylation_50kb.freq"))
#allcov <- getCoverage(bis)

#fifth <- quantile(allcov, probs = 0.05)
#ninety_fifth <- quantile(allcov, probs = 0.95)
# filter for for calls between the 5th and 95th perentile for coverage

#meth <- getMeth(bis, type = "raw")
#getmeth <- dplyr::as_tibble(cbind(freq$X1, freq$X2, meth,allcov))%>%
#  dplyr::rename("chrom" = 1, "start" = 2, "meth" = 3, "cov" = 4) %>%
#  dplyr::mutate("cov" = as.numeric(cov)) %>%
# # filter(cov < ninety_fifth) %>%
# # filter(cov >  fifth) %>%
#  mutate("meth" = as.numeric(meth)) %>%
#  mutate("start" = as.numeric(start)) %>%
#  mutate("chrom" = as.factor(chrom)) 


freq.gr <- GRanges(freq)
saveRDS(freq.gr, file =paste0(dat, "/methylation_calls/chm13_methylation_NanopolishFreq_50kb.rds"))

keepi <- findOverlaps(freq.gr,chm13)
freq.matched <- freq.gr[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
    mcols(freq.matched),
    mcols(chm13[subjectHits(keepi)]))

freq.matched

final_dat <- as.data.frame(freq.matched) %>%
  distinct() %>%
  GRanges()
#%>%
#  group_by(seqnames, start, end, meth) %>%
#  summarise(CpG = mean(CpG))

saveRDS(final_dat, file =paste0(dat, "/methylation_calls/chm13_methylation_200bp_slidingCpG.rds"))
