library(tidyverse)
library(zoo)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

HG01109 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG01109_methylation_frequency.tsv") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(cell_line = "HG01109") %>%
  mutate(norm_start = start - min(start))

HG01243 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG01243_methylation_frequency.tsv") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(cell_line = "HG01243")%>%
  mutate(norm_start = start - min(start))

HG03098 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG03098_methylation_frequency.tsv") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(cell_line = "HG03098")%>%
  mutate(norm_start = start - min(start))

HG03492 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG03492_methylation_frequency.tsv") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(cell_line = "HG03492") %>%
  mutate(norm_start = start - min(start))

HG002 <- read_tsv(paste0(dat,"/revision_analysis/HG002_pooled/HG002_HPRC_pooledCpGmethylationFrequency_HG002chrX.tsv")) %>%
  filter(chromosome=='chrX') %>%
  filter(start > 55941025) %>%
  filter(end < 59053403) %>%
  mutate(cell_line = "HG002") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(norm_start = start - min(start))

CHM13 <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  filter(chromosome=='chrX') %>%
  filter(start> 57820108 ) %>%
  filter(end < 60927026) %>%
  mutate(cell_line = "CHM13") %>%
  mutate(smoothed = rollmean(methylated_frequency, 300, fill = NA)) %>%
  mutate(norm_start = start - min(start))

meth_dat <- rbind(HG01109,HG01243,HG03098,HG03492, HG002, CHM13) 
sub <- meth_dat[seq(1, nrow(meth_dat), 100), ]
ggplot(data=sub, aes(x=norm_start/1e6,y=smoothed))+geom_line(stat="identity")+facet_wrap(~cell_line, ncol=1)+theme_classic()+labs(y="Methylation Frequency", x="[Mb]")+ylim(0,1)

ggsave(
  paste0(figs, "/methyl_profiles/PanXMethylation.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 8, 
  useDingbats=FALSE
)

#### k2 annotation panes

k2.hg01109 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG01109_cenX_hap_k2.nc.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "HG01109") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

k2.hg0002 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG002_june_cenX_hap_k2.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "HG002") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

k2.hg01243 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG01243_cenX_hap_k2.nc.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "HG01243") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

k2.hg03098 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG03098_cenX_hap_k2.nc.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "HG03098") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

k2.hg03492 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/HG03492_cenX_hap_k2.nc.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "HG03492") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

k2.chm13 <-  read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HPRC_CDRs/CHM13_june_cenX_hap_k2.bed", col_names = c("chr", "start", "end", "hap")) %>%
  mutate(cell_line = "CHM13") %>%
  mutate(norm_start = start - min(start), len=end-start, norm_end=norm_start+len)

all.haps <-rbind(k2.hg01109,k2.hg0002,k2.hg01243,k2.hg03098,k2.chm13,k2.hg03492)
ggplot(data=all.haps, aes(xmin=norm_start/1e6,xmax=norm_end/1e6, ymin=0,ymax=1, color=as.factor(hap),fill=as.factor(hap)))+geom_rect()+facet_wrap(~cell_line, ncol=1)+theme_classic()+labs(y="Methylation Frequency", x="[Mb]")+facet_wrap(~cell_line, ncol=1)+theme_classic()+labs(y="Methylation Frequency", x="[Mb]")+ylim(0,1)+scale_color_manual(values = c("1" = "darkred", "2" = "gray"))+theme(legend.position = "none")

ggsave(
  paste0(figs, "/methyl_profiles/PanXHap2.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 3, 
  useDingbats=FALSE
)
