#!/usr/bin/Rscript

library(bsseq)
library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(BiocParallel)

path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/pooled"
files=list.files(path, pattern= ".bismark.out")


bisobj <- read.bismark(paste0(path,"/",files))

saveRDS(bisobj, paste0(path, "/HG002chrX_nanonome_GpC_merged_unsmoothedbsojb.rds"))                     
#bsobj.merged <- collapseBSseq(bisobj,group = c("merged", "merged", "merged"))

bsobj.smoothed <-  BSmooth(bisobj,h=100,ns=10,maxGap=10^3,
                           BPPARAM = MulticoreParam(workers = 10,progressbar = TRUE),
                           verbose=TRUE)

saveRDS(bsobj.smoothed, paste0(path, "/HG002chrX_nanonome_GpC_Smoothedbsojb.rds"))       


Cov <- getCoverage(bsobj.smoothed)
keepLoci = (Cov >= 5)
Smoothed.subset = bsobj.smoothed[keepLoci, ]

minwin = 50
a = 0.01
cutoff = NULL
qcutoff = 0.99

peaks <- gpcPeakCaller(Smoothed.subset)

write_tsv(peaks, paste0(path, "/HG002chrX_Peaks.tsv"),quote=F)
