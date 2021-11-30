library(tidyverse)
library(ggpubr)
library(ggridges)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

hg002_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
  # mutate(called_sites_unmethylated=called_sites-called_sites_methylated) %>%
  GRanges()

res <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/walusat/WAluSat_chm13v1.0_mon-array.txt") %>%
  mutate(ID=row_number()) %>%
  GRanges()

ovls <- FindOvls(res,chm13_meth)

ovls.stat <- ovls %>%
  group_by(ID, group,seqnames) %>%
  summarise(methylation=mean(methylated_frequency))

my_comparisons <- list(c("ARRAY", "MONOMER"), c("ARRAY", "DUP"),c("MONOMER","DUP"))


p1 <- ovls.stat %>%
 # unite("name", single_or_array, centromeric) %>%
 # filter(name != "Singleton_nonCEN") %>%
ggplot(aes(x=group,y=methylation, fill=group))+geom_violin(width=1)+geom_boxplot(width=.1, outlier.shape = NA)+ stat_compare_means(comparisons = my_comparisons, test = "kruskal.test")+theme_classic()

ggsave(filename = paste0(figs, "/CHM13_WALUSAT_MethylationArrayvsSingleton.pdf"), 
       plot = p1,
       width = 8,
       height = 5)




ggsave(filename = paste0(figs, "/CHM13_LSAU_MethylationBSatvsD4Z4.pdf"), 
       plot = p2,
       width = 8,
       height = 5)

ovls <- FindOvls(res,hg002_meth)

ovls.stat <- ovls %>%
  filter(called_sites > 10) %>%
  group_by(ID, group,seqnames,start,end) %>%
  summarise(methylation=mean(methylated_frequency))

p3<- ovls.stat %>%
  #unite("name", single_or_array, centromeric) %>%
  #filter(name != "Singleton_nonCEN") %>%
  ggplot(aes(x=group,y=methylation, fill=group))+geom_violin(width=1)+geom_boxplot(width=.1, outlier.shape = NA)+ stat_compare_means(comparisons = my_comparisons, method = "t.test")+theme_classic()

ggsave(filename = paste0(figs, "/HG002_LSAU_MethylationArrayvsSingleton.pdf"), 
       plot = p3,
       width = 8,
       height = 5)

p4 <- ovls.stat %>%
  unite("name", single_or_array, centromeric,seqnames) %>%
  filter(name %in% c("array_nonCEN_chr4", "array_CEN_chr1")) %>%
  ggplot(aes(x = name, y = methylation, fill=name)) +
  geom_dotplot(binaxis = "y", stackdir = "center")+ stat_compare_means( method = "t.test")+theme_classic()


ggsave(filename = paste0(figs, "/HG002_LSAU_MethylationBSatvsD4Z4.pdf"), 
       plot = p4,
       width = 8,
       height = 5)
