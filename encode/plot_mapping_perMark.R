

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)

# load data
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

chm13_en <- read_csv(paste0(dat, "/encode/210308_alignments/mapping_stats_reordered.txt"), skip=1) %>%
  select(c(1:8)) %>%
  dplyr::rename("Total_reads" = `Total Reads`, "Cell_type"= `Cell type`) %>%
  gather(mapping, chm13,-c(Cell_type,Experiment,Target)) 

hg38_en <- read_csv(paste0(dat, "/encode/210308_alignments/mapping_stats_reordered.txt"), skip=1) %>%
  select(c(1:4,9:12)) %>%
  dplyr::rename("Total_reads" = `Total Reads`, "Cell_type"= `Cell type`, "mapped"=mapped_1, "filtered"=filtered_1,"dedup"=dedup_1,"kmer"=kmer_1) %>%
  gather(mapping, hg38,-c(Cell_type,Experiment,Target)) 

en=merge(chm13_en,hg38_en,by=c("Experiment","Cell_type","Target","mapping")) %>%
  filter(mapping != "Total_reads") %>%
  group_by(mapping,Target) %>%
  summarize(perc=((chm13/hg38)-1)*100) 


en$mapping <- factor(en$mapping, levels = c("mapped", "filtered", "dedup", "kmer"))
ggplot(en, aes(y=perc,x=mapping, color=Target))+geom_boxplot()+theme_classic()+geom_dotplot(aes(fill=Target),binaxis='y', stackdir='center', dotsize=.2,position=position_dodge(.75))


en.stat= merge(chm13_en,hg38_en,by=c("Experiment","Cell_type","Target","mapping")) %>%
  filter(mapping != "Total_reads") %>%
  group_by(mapping,Target,Cell_type) %>%
  summarize(perc=((chm13/hg38)-1)*100) %>%
  group_by(mapping,Target,Cell_type) %>%
  summarise(sd=sd(perc),med=median(perc))

en.stat$mapping <- factor(en.stat$mapping, levels = c("mapped", "filtered", "dedup", "kmer"))
ggplot(en.stat, aes(y=med,x=mapping, color=Target,group=Target))+geom_line()+theme_classic()+geom_point(size=2.5)+facet_wrap(~Cell_type)

ggsave(
  paste0(figs, "/", "ENCODE_FacetedCellType_lineplot.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 10,
)

lines=c("HAP-1", "SJSA1", "SJCRH30")


en.stat= merge(chm13_en,hg38_en,by=c("Experiment","Cell_type","Target","mapping")) %>%
  filter(mapping != "Total_reads") %>%
  filter(Cell_type %in% lines) %>%
  group_by(mapping,Target,Cell_type) %>%
  summarize(perc=((chm13/hg38)-1)*100) %>%
  group_by(mapping,Target,Cell_type) %>%
  summarise(sd=sd(perc),med=median(perc))

en.stat$mapping <- factor(en.stat$mapping, levels = c("mapped", "filtered", "dedup", "kmer"))
ggplot(en.stat, aes(y=med,x=mapping, color=Target,group=Target))+geom_line()+theme_classic()+geom_point(size=2.5)+facet_wrap(~Cell_type,ncol=1)+ylim(0,3)


en.stat= merge(chm13_en,hg38_en,by=c("Experiment","Cell_type","Target","mapping")) %>%
  filter(mapping != "Total_reads") %>%
  group_by(mapping,Target) %>%
  summarize(perc=((chm13/hg38)-1)*100) %>%
  group_by(mapping,Target) %>%
  summarise(sd=sd(perc),med=median(perc))

en.stat$mapping <- factor(en.stat$mapping, levels = c("mapped", "filtered", "dedup", "kmer"))
ggplot(en.stat, aes(y=med,x=mapping, group=Target))+geom_line(aes(color=Target))+theme_classic()+theme(text = element_text(size=20))+labs(y="Percent Change in CHM13", x="Mapping Step")+geom_point(aes(color=Target),size=3,shape="square")

ggsave(
  paste0(figs, "/", "ENCODE_meta_lineplot.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)
