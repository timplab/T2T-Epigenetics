source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(BSgenome.HG002.chrX)
library(zoo)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002"

hg002_meth <- read_csv(paste0(dat, "/nanonome/methylation_calls/pooled/HG002_nanonome_chrX_CpGmethylation_ByCall.csv")) %>%
  GRanges() 

cenpb <- read_delim(paste0(dat,"/reference/chrx_fuzznuc.bed"), delim=" ") %>%
  rename(`  Start` = "start",`    End`="end", ` Strand` = "strand" ) %>%
  select(c("start", "end", "strand")) %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  filter(start > 54000000) %>%
  filter(end < 62000000) %>%
  mutate(reg=case_when(start > 57490000 & end < 57570000 ~ "CDR", 
                       TRUE ~"non-CDR")) %>%
  mutate(seqnames="chrX") %>%
  mutate(strand=case_when(strand=="      +" ~ "+", 
                          TRUE ~ "-")) %>%
  GRanges()
     
cenpb.meth <- FindOvls(cenpb, hg002_meth) %>%
  group_by(reg,start) %>%
  summarize(methylated_frequency=mean(mcall))

ggplot(data=cenpb.meth, aes(x=reg, y=methylated_frequency))+geom_boxplot()

kruskal.test(methylated_frequency ~ reg, data = cenpb.meth)

cdr <- cenpb.meth %>%
  filter(reg == "CDR")
non_cdr <- cenpb.meth %>%
  filter(reg != "CDR")

ggsave(
  paste0(figs, "/CDR_vsnonCDR_cenpb.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 4,
  height = 5
)
