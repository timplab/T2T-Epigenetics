#!/usr/bin/Rscript

knitr::opts_chunk$set(echo = FALSE, warning=FALSE)
library(knitr)
library(tidyverse)
library(bsseq)
library(Biostrings)
library(ggplot2)
library(png)
library(cowplot)
options(scipen=999)
library(zoo)
library(ggpmisc)
library(Repitools)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls"


getQQ <- function(one,two,qs = seq(0,1,0.01)) {
  tibble(one = quantile(one,qs,na.rm = T),
         two = quantile(two,qs,na.rm = T),
         qunatile = qs)
}
# https://stackoverflow.com/questions/43875716/find-start-and-end-positions-indices-of-runs-consecutive-values
#getRuns <- function(calls){
#  calls <- calls %>%
#    filter(mcall != -1) %>%
#    group_by(qname)
#  calls.list <- calls %>%
#    group_split(keep = F)
#  names(calls.list) <- calls %>% group_keys() %>% .$qname
#  runs.list <- lapply(calls.list,function(x){
#    rle(x$mcall) %>%
#    unclass() %>% as.tibble() %>%
#    mutate( endi = cumsum(lengths),
#            starti = c(1,dplyr::lag(endi)[-1]+1),
#            start = x$start[starti],
#            end = x$start[endi],
#            width = end - start + 1) %>%
#      filter( width >= 0) # remove negative widths (in case of dups, etc.)
#  })
#  runs <- bind_rows(runs.list,.id = "qname")
#}
order_reads <- function(x,bounds=NULL){
  # get boundaries of reads if not provided
  if (is.null(bounds)){
    bounds <- x%>% group_by(qname) %>%
      summarize(start = min(start),
                end = max(end))
    # label y based on order of smallest start
    bounds<- bounds %>% arrange(start, end) %>%
      mutate(
        readi = seq_len(length(unique(x$qname))),
        ymin = -readi - 0.8, 
        ymax = ymin + 0.6)
  }
  x <- x %>%
    mutate(ymin = bounds$ymin[match(qname,bounds$qname)],
           ymax = bounds$ymax[match(qname,bounds$qname)])
  bounds <- bounds %>%
    mutate(ymin = bounds$ymin[match(qname,bounds$qname)],
           ymax = bounds$ymax[match(qname,bounds$qname)])
  return(list(x = x,bounds = bounds))
}

smoothCalls <- function(calls,reg=NULL,bandwidth = 80){
  calls <- calls %>%
    mutate(mcall = ifelse(abs(score)>1,sign(score),score)) # ceiling based on log-lik ratio - this uses log-lik ratio when call is not obvious
  if (is.null(reg)) {
    xpoints <- seq(min(calls$start),max(calls$start))
  } else {
    reg <- as_tibble(reg)
    xpoints <- seq(reg$start,reg$end)
  }
  ks <- ksmooth(calls$start,calls$mcall,bandwidth = 80,kernel = "normal",x.points = xpoints)
  tibble(
    start = ks$x,
    mcall_smooth = ks$y, 
    mcall = case_when(
      mcall_smooth > 0 ~ 1,
      mcall_smooth < 0 ~ 0,
      TRUE ~ -1)) 
}

readsGC <- tabix_mbed(paste0(dat, "/HG002_nanonome_chrX_GpCmethylation.merge1-3.bed.gz"),extcol = "motif",by = "read") 

size_sel <- readsGC %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 20000) %>%
  mutate(motif = "GC")

gccalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)


calls.reg <- gccalls %>%
 # filter(start > 55000000) %>%
 # filter(end < 60000000) %>%
    group_by(qname)
    
group_names <- group_keys(calls.reg)$qname

calls.list <- calls.reg %>% 
group_split(keep = T)

names(calls.list) <- group_names
smooth.list <- lapply(calls.list,smoothCalls)
calls.smooth <- bind_rows(smooth.list,.id = "qname")
cgsmooth.list <- calls.smooth

calls.gr <-  calls.smooth %>%
  ungroup() %>%
  group_by(start) %>%
  summarise(num_meth = sum(mcall == 1), num_unmeth = sum(mcall == 0)) %>%
  mutate(chr = "chrX", start = start, end = start) %>%
 # filter(start > 112250000) %>%
 # filter(end < 112550000) %>%
  GRanges()

saveRDS(calls.gr, file=paste0(dat, "/GpC_smoothed_calls_20kb_80bp.rds"))



