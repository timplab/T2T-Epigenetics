library(tidyverse)
library(ggpubr)
library(ggridges)
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated)

hg002_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv") %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated)


array <- read_csv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/LSAU-BSAT-composite_chm13v1.0_loci.csv", skip=1, col_names=c("chr", "start", "end", "name", "centromeric", "len", "single_or_array")) %>%
  mutate(ID=row_number()) %>%
  filter(single_or_array=="array") %>%
  GRanges()

LSAU <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_RepeatMaskerV2_LSAU-BSAT.bed", col_names=c("chrom",  "start",      "end",        "name",    "score",   "strand",  "thickStart",      "thickEnd",        "reserved" ,       "swScore", "repClass" ,       "repSubClass" ,    "repDiverence")) %>%
  GRanges()


lsau_arrays <- FindOvls(LSAU,array) %>%
  select(c(seqnames, start, end, name, centromeric, len, single_or_array)) 

singles <- read_csv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/LSAU-BSAT-composite_chm13v1.0_loci.csv", skip=1, col_names=c("seqnames", "start", "end", "name", "centromeric", "len", "single_or_array")) %>%
  filter(single_or_array=="Singleton")

all.lsau <- rbind(singles, lsau_arrays) %>%
  mutate(ID=row_number()) %>%
  GRanges()

ovls <- FindOvls(GRanges(chm13_meth),all.lsau)

ovls.stat <- ovls %>%
  group_by(ID, centromeric, single_or_array,seqnames) %>%
  summarise(methylation=mean(methylated_frequency))

my_comparisons <- list(c("array_CEN", "Singleton_CEN"), c("array_nonCEN", "Singleton_CEN"))


p1 <- ovls.stat %>%
  unite("name", single_or_array, centromeric) %>%
  filter(name != "Singleton_nonCEN") %>%
ggplot(aes(x=name,y=methylation, fill=name))+geom_violin(width=1)+geom_boxplot(width=.1, outlier.shape = NA)+ stat_compare_means(comparisons = my_comparisons, method = "t.test")+theme_classic()

ggsave(filename = paste0(figs, "/CHM13_LSAU_MethylationArrayvsSingleton.pdf"), 
       plot = p1,
       width = 8,
       height = 5)


p2 <- ovls.stat %>%
  unite("name", single_or_array, centromeric,seqnames) %>%
  filter(name %in% c("array_nonCEN_chr4", "array_CEN_chr1")) %>%
  ggplot(aes(x = name, y = methylation, fill=name)) +
  geom_dotplot(binaxis = "y", stackdir = "center")+ stat_compare_means( method = "t.test")+theme_classic()


ggsave(filename = paste0(figs, "/CHM13_LSAU_MethylationBSatvsD4Z4.pdf"), 
       plot = p2,
       width = 8,
       height = 5)

ovls <- FindOvls(all.lsau,GRanges(hg002_meth))

ovls.stat <- ovls %>%
  filter(called_sites > 10) %>%
  group_by(ID, centromeric, single_or_array,seqnames,start,end) %>%
  summarise(methylation=mean(methylated_frequency))

p3<- ovls.stat %>%
  unite("name", single_or_array, centromeric) %>%
  filter(name != "Singleton_nonCEN") %>%
  ggplot(aes(x=name,y=methylation, fill=name))+geom_violin(width=1)+geom_boxplot(width=.1, outlier.shape = NA)+ stat_compare_means(comparisons = my_comparisons, method = "t.test")+theme_classic()

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
