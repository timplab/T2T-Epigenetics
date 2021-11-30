#!/usr/bin/Rscript

library(bsseq)
library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(BiocParallel)

path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome"
files=list.files(path, pattern= "_bismark.out")


rows=c("run1", "run2", "run3")
columndat = data.frame(rows,paste0(path,"/",files)) %>%
  column_to_rownames("rows")

bisobj <- read.bismark(paste0(path,"/",files),
                       colData=columndat)

saveRDS(bisobj, paste0(path, "/HG002_nanonome_GpC_unmerged_unsmoothedbsojb.rds"))                     
bsobj.merged <- collapseBSseq(bisobj,group = c("merged", "merged", "merged"))

bsobj.smoothed <-  BSmooth(bsobj.merged,h=100,ns=10,maxGap=10^3,
                           BPPARAM = MulticoreParam(workers = 10,progressbar = TRUE),
                           verbose=TRUE)

saveRDS(bsobj.smoothed, paste0(path, "/HG002_nanonome_GpC_Smoothedbsojb.rds"))       


Cov <- getCoverage(bsobj.smoothed)
keepLoci = (rowSums(Cov >= 5) >= 3)
Smoothed.subset = bsobj.smoothed[keepLoci, ]

minwin = 50
a = 0.01
cutoff = NULL
qcutoff = 0.99

gpcPeakCaller <- function(bsobj,minwin = 50,a = 0.01,cutoff = NULL, qcutoff = 0.99, samp=1){
  gpc.meth <- getMeth(bsobj,type="smooth",what="perBase")[,samp]
  # remove NA
  keepi <- !is.na(gpc.meth)
  bsobj <- bsobj[keepi,]
  gpc.meth <- gpc.meth[keepi]
  gpc.cov <- getCoverage(bsobj,type="Cov",what="perBase")[,samp]
  gpc.m <- getCoverage(bsobj,type="M",what="perBase")[,samp]
  ## compare to baseline ----
  baseline <- median(gpc.meth,na.rm = T)
  gpc.diff <- gpc.meth - baseline
  if ( is.null(cutoff)) {
    cutoff <- quantile(gpc.diff,qcutoff,na.rm = T)
  }
  gpc.direction <- ifelse(gpc.diff > cutoff, 1, -1) # cut off by qcutoff
  gpc.gr <- granges(bsobj)
  chrs <- as.character(seqnames(gpc.gr))
  pos <- start(gpc.gr)
  ## find regions of + ----
  regions <- bsseq:::regionFinder3(gpc.direction,chrs,pos)$up
  regions <- as_tibble(regions)
  ## then add average and peak accesibility, along with coverage, then binomial test ---
  regions <- regions %>%
    rowwise() %>%
    mutate(
      coverage = sum(gpc.cov[idxStart:idxEnd]),
      methylated = sum(gpc.m[idxStart:idxEnd]),
      average = mean(gpc.meth[idxStart:idxEnd]),
      peak = max(gpc.meth[idxStart:idxEnd]),
      p.value = binom.test(methylated,coverage,baseline,alternative = "greater")$p.value
    ) %>%
    ungroup()
  ## multiple test adjusting using FDR
  regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")
  
  ## significance based on :
  ## 1. width
  ## 2. alpha
  ## 3. minimum peak height
  minpeak <- 0 #baseline + 2 * cutoff 
  regions %>%
    mutate(
      width = end - start + 1,
      sig = ifelse(
        adjusted.pval <= a &
          width >= minwin &
          peak >= minpeak
        ,"sig","insig"))
}

peaks1 <- gpcPeakCaller(Smoothed.subset,samp=1)
peaks2 <- gpcPeakCaller(Smoothed.subset,samp=2)
peaks3 <- gpcPeakCaller(Smoothed.subset,samp=3)


write_tsv(peaks1, paste0(path, "/HG002_nanonome1_winnowmapk15_Peaks.tsv"),quote=F)
write_tsv(peaks2, paste0(path, "/HG002_nanonome2_winnowmapk15__Peaks.tsv"),quote=F)
write_tsv(peaks3, paste0(path, "/HG002_nanonome4-5_winnowmapk15__Peaks.tsv"),quote=F)
