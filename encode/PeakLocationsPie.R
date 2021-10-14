library(tidyverse)
library(GenomicRanges)
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")


FindOvls <- function(obj1,obj2){
  
  keepi <- findOverlaps(obj1,obj2, type = "any")
  freq.matched <- obj1[queryHits(keepi)]
  
  mcols(freq.matched) <- cbind.data.frame(
    mcols(freq.matched),
    mcols(obj2[subjectHits(keepi)]))
  return(as.data.frame(freq.matched))
}

bed <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/pie/repeat_all.bed", col_names=c("chr", "start", "end", "name", "region")) %>%
  GRanges()

genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/pie/genes_all.bed", col_names=c("chr", "start", "end", "name", "score", "strand", "geneStart", "geneEnd", "score2", "num1", "num2", "nums", "gene_id", "none", "col1", "col2", "col3", "col4", "col5", "region")) %>%
  dplyr::select(c(chr, start,end,name,region)) %>%
  GRanges()

marks <- c("H3K36me3", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me3","CTCF")

df <- data.frame()

for (mark in marks){
  peaks <- read_tsv(paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks/cat/", mark, "_unique.bed"),col_names = c("chr", "start", "end", "samp", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
    mutate(mark = mark) 
  df <- rbind(df, peaks)
  
}

df.gr <- GRanges(df)

nonovl <- as.data.frame(df.gr[!df.gr %over% bed,]) %>%
  mutate(name = "none", region = "other") %>%
  dplyr::select(-c(name)) %>%
  distinct()
  


ovl <- FindOvls(df.gr, bed) %>%
  as.data.frame()%>%
  dplyr::select(-c(name)) %>%
  distinct()

all.dat <- rbind(nonovl,ovl) %>%
  distinct() %>%
  group_by(mark) %>%
  mutate(total=n()) %>%
  group_by(region, mark) %>%
  summarise(n=n(), percent = (n/total)*100, total=total) %>%
  distinct()

ggplot(data=all.dat, aes(x=1,percent, fill = region))+geom_bar(stat="identity", color="white")+
  coord_polar("y", start=0) +
  theme_void()+facet_wrap(~mark)

ggsave('/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/pie/Pie.pdf', 
       last_plot(), 
       height = 8, 
       width=10)


ovl <- FindOvls(df.gr, genes) %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, name, mark,region) %>%
  distinct()

all.dat <- ovl %>%
  group_by(mark) %>%
  mutate(total=n()) %>%
  group_by(region, mark) %>%
  summarise(n=n(), percent = (n/total)*100, total=total) %>%
  distinct()
  
ggplot(data=all.dat, aes(x=1,percent, fill = region))+geom_bar(stat="identity", color="white")+
  coord_polar("y", start=0) +
  theme_void()+facet_wrap(~mark)

ggsave('/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/pie/GenesPie.pdf', 
       last_plot(), 
       height = 8, 
       width=10)
