source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(GenomicRanges)
library(Repitools)
library(rtracklayer)

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

SAT = c("BSAT","HSAT1", "HSAT2", "HSAT3","HOR", "MON", "CT")

censat = read_tsv(paste0(dat, "/revision_analysis/annotations/pie/chm13_v1_MarkerDesertsIntersectRepeat_all.bed"), col_names = c("chr", "start", "end", "name", "region")) %>%
#  mutate(name = ifelse(grepl("gap", name), "gap-rDNA", name)) %>%
#  mutate(name = ifelse(grepl("bsat", name), "BSat", name)) %>%
#  mutate(name = ifelse(grepl("hsat1", name), "HSat1", name)) %>%
#  mutate(name = ifelse(grepl("hsat2", name), "HSat2", name)) %>%
#  mutate(name = ifelse(grepl("hsat3", name), "HSat3", name)) %>%
#  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
#  mutate(name = ifelse(grepl("mon", name), "MON", name)) %>%
#  mutate(name = ifelse(grepl("ct", name), "CT", name)) %>%
#  mutate(newname = ifelse(region == "sd_only", name, "sd")) %>%
#  mutate(newname = ifelse(region == "cen_sd", name, "sd")) %>%
  dplyr::select(-c(name)) %>%
  distinct() 

stat.sum <- censat %>%
  mutate(len=end-start) %>%
  group_by(region) %>%
 # dplyr::filter(name %in% SAT) %>%
  summarise(avg = mean(len), sum = sum(len))

ggplot(stat.sum, aes(x="",y=(sum/3056916522)*100, fill=region))+
  geom_bar(stat="identity", color="black")+labs(x="", y="Mb")+
  theme(text = element_text(size=18))+
  coord_polar("y", start=0) +
  theme_void()+
  geom_text(aes(label = (round((sum/3056916522)*100,2))),
            position = position_stack(vjust = 0.5))+
  scale_fill_brewer(palette="Set1")

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures/KmerDesertPie.pdf",
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 5,
)
