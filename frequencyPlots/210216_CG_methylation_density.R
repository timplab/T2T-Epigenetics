#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(bsseq)
library(Biostrings)
library(ggplot2)
library(png)
library(cowplot)
options(scipen=999)
library(zoo)
library(BSgenome.t2t.v1.0.release)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

# set data output path

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/violin"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load and parse censat files

censat = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("hsat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("hsat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("hsat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("hsat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("hsat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("hsat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("dhor", X4), "DHOR", name)) %>%
  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
  dplyr::select(c(X1, X2, X3, name)) %>%
  dplyr::rename("chrom" =1, "start" = 2 ,"end" =3) 

#write.table(censat, paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), sep="\t", col.names = F, row.names = F, quote=F)
# set repeat colors
censatColors =c("(CATTC)n" = "#E87C71",
                "(GAATC)n"="#E28455",
                "HOR"="#D78C32",
                "BSAT"="#E370AB",
                "CER" = "#CE9334",
                "HSAT2"="#C19935",
                "HSAT1"="#A2A638",
                "HSAT3"="#8CAC3E",
                "Low_complexity"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "MST"="#53BB73",
                "GSAT"="#55BE8D",
                "RNA"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "ACRO1"="#9400D3",
                "HSAT4"="#53B0E3",
                "SATR"="#5AA5DA",
                "CT"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "MON"="#A885BC",
                "SST"="#C378B2",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "gap-rDNA"="#ff4000",
                "TE" = "#ffbf00", 
                "TAR"= "#0080ff",
                "ACRO"="#9400D3", 
                "DHOR" = "gray")

# load CG and methylation GRanges data 
chm13_CpG <- readRDS(paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))
#chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_NanopolishFreq_50kb.rds"))

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()


# find overlaps between censat regions and methylation calls
censat.gr <- GRanges(censat)
keepi <- findOverlaps(chm13_meth,censat.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

# split into groups for plotting

SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT")



censat_meth <- as.data.frame(freq.matched) %>%
filter(name %in% SAT)

# violing plot of methylation frequency

stats <- censat_meth %>%
  group_by(name) %>%
  summarize(med = median(methylated_frequency))

violin <- ggplot(data = censat_meth, aes(x = factor(name), y = methylated_frequency, fill = name))+geom_violin(width=1.5)+theme_classic(base_size = 12)+ scale_fill_manual(values = censatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")+geom_boxplot(width=0.1,outlier.shape = NA)+geom_hline(yintercept = .3680, linetype= "dashed") 

ggsave(
  paste0(figs, "/", "AllSatv2.021621_methylation_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)

for (i in SAT){
  print(i)
  meth <- censat_meth %>%
    filter(name == i)
violin <- ggplot(data = meth, aes(y = methylated_frequency, fill = name, x=seqnames))+geom_violin()+theme_classic(base_size = 20)+ scale_fill_manual(values = censatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")+geom_boxplot(width=0.1,outlier.shape = NA)+geom_hline(yintercept = .3680, linetype= "dashed")

ggsave(
  paste0(figs, "/", i, "_methylationv2.021621_freqByChr.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)

}


keepi <- findOverlaps(chm13_CpG,censat.gr)
freq.matched <- chm13_CpG[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))



cpg <- as.data.frame(freq.matched) %>%
    filter(name %in% SAT)
  
stats.cpg <- cpg %>%
  group_by(name) %>%
  summarize(med = median(CpG))

  # violing plot of methylation frequency

violin <- ggplot(cpg, aes(x = factor(name), y = CpG, fill = name))+ scale_fill_manual(values = censatColors, drop = FALSE)+theme_classic(base_size = 20)+labs(x = "Repeat", y = "CpG Density")+geom_boxplot(outlier.shape = NA)+coord_cartesian(y=c(0,20))+geom_hline(yintercept = 3, linetype= "dashed") 
  
  ggsave(
    paste0(figs, "/AllSatv2.021621_CpG_boxplot.pdf"),
    plot = violin,
    scale = 1,
    width = 8,
    height = 5,
  )
  
  for (i in SAT){
    print(i)
    cpg.sub <- cpg %>%
      filter(name == i)
    
    violin <- ggplot(cpg.sub, aes(x = seqnames, y = CpG, fill = name))+ scale_fill_manual(values = censatColors, drop = FALSE)+theme_classic(base_size = 20)+labs(x = "Repeat", y = "CpG Density")+geom_boxplot(outlier.shape = NA)+coord_cartesian(y=c(0,20))+geom_hline(yintercept = 3, linetype= "dashed")
    
    ggsave(
      paste0(figs, "/", i, "_CpG_v2.021621freqByChr.pdf"),
      plot = violin,
      scale = 1,
      width = 8,
      height = 5,
    )
    
  }
  
 