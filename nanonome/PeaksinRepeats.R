library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(ggbreak)
options(scipen = 100)

list=c("HSAT2", "HSAT3", "HSAT4", "BSAT", 'GSAT', "ACRO","SST1_Composite", "MON", "DHOR", "HSAT1")


censat.peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/chm13_v1_CenSat_PeakIntersections.bed", col_names=c("chr", "peak_start", "peak_end", "num", "strand1", "idx1", "idx2", "cov", "num1", "num2","num3", "meth", "meth2", "p", "padj", "sig", "rep_chr", "rep_start", "rep_end", "name"))

peaks.gr <- censat.peaks %>%
  select(c(chr,peak_start,peak_end)) %>%
  GRanges()

peaks.gr <- peaks.gr[!peaks.gr %over% hg002.flag,]


censat.widths <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
  mutate(name = ifelse(grepl("hsat1", name), "HSAT1", name)) %>%
  mutate(name = ifelse(grepl("hsat2", name), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("hsat3", name), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("hsat4", name), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("hsat5", name), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("hsat6", name), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", name), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", name), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("dhor", name), "DHOR", name)) %>%
  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", name), "MON", name)) %>%
  mutate(name = ifelse(grepl("GSAT", name), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("ACRO", name), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("SST1_Composite", name), "SST1_Composite", name)) %>%
  mutate(width=chromEnd-chromStart) %>%
  filter(name %in% list) %>%
  group_by(name) %>%
  summarise(width=sum(width))


censat.peaks.stat <- censat.peaks %>%
  mutate(name = ifelse(grepl("hsat1", name), "HSAT1", name)) %>%
  mutate(name = ifelse(grepl("hsat2", name), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("hsat3", name), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("hsat4", name), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("hsat5", name), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("hsat6", name), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", name), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", name), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("dhor", name), "DHOR", name)) %>%
  mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", name), "MON", name)) %>%
  mutate(name = ifelse(grepl("GSAT", name), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("ACRO", name), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("SST1_Composite", name), "SST1_Composite", name)) %>%
  mutate(width=rep_end-rep_start) %>%
  filter(name %in% list) %>%
  group_by(name) %>%
  summarize(n=n()) 


censat.merge <- merge(censat.peaks.stat, censat.widths) %>%
  mutate(peaks=(n/width)*1e6)

rep.peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/chm13_v1_RepeatMaskerV2_PeakIntersections.bed", col_names=c("chr", "peak_start", "peak_end", "num", "strand1", "idx1", "idx2", "cov", "num1", "num2","num3", "meth", "meth2", "p", "padj", "sig", "rep_chr", "rep_start", "rep_end", "name", "rep_score","rep_strand",  "num4", "num5", "nums", "len", "repClass", "rep_class"))


list2 <- c("LINE", "SINE", 'LTR', "Satellite")
rep.peaks.stat  <- rep.peaks %>%
  mutate(width=rep_end-rep_start) %>%
  group_by(repClass) %>%
  summarize(n=n()) 

rep.widths  <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_RepeatMaskerV2.bed") %>%
  mutate(width=chromEnd-chromStart) %>%
  filter(repClass %in% list2) %>%
  group_by(repClass) %>%
  summarize(width=sum(width)) 

rep.merge <- merge(rep.peaks.stat, rep.widths) %>%
  mutate(peaks=(n/width)*1e6) %>%
  select(c(repClass,n,peaks))


wg <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome/Peaks/Cat_Peaks_BlacklistFiltered.bed", col_names=c("chr", "start", "end")) %>%
  mutate(repClass="WG") %>%
  group_by(repClass) %>%
  summarise(n=n())  %>%
  mutate(peaks=(n/3056916522)*1e6)

rep.merge <- rbind(rep.merge,wg)

ggplot(rep.merge, aes(x= reorder(repClass, peaks), y=peaks))+geom_bar(stat="identity")+coord_flip()
ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/AllRepeatsBar.pdf", 
       last_plot(), 
       height = 5,
       width=5)

ggplot(censat.merge, aes(x=reorder(name, peaks), y=peaks))+geom_bar(stat="identity")+coord_flip()

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/SatRepeatsBar.pdf", 
       last_plot(), 
       height = 5,
       width=5)

sst1.peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/SST1_CEN-array-mon_NONCEN-array_withLabels_PeakIntersections.bed", col_names=c("chr", "peak_start", "peak_end", "rep_chr", "rep_start", "rep_end", "rep_strand", "group"))%>%
  group_by(group) %>%
  summarize(n=n()) 

sst1.all <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/revision/SST1_CEN-array-mon_NONCEN-array_withLabels.txt")%>%
  mutate(width=end-start) %>%
  group_by(group) %>%
  summarise(width=sum(width))

sst1.merge  <- merge(sst1.peaks,sst1.all) %>%
  filter(group != "mon_CEN") %>%
  mutate(peaks=(n/width)*1e6)

censat.merge <- censat.merge %>%
  rename("group"=name) %>%
  filter(group != "SST1_Composite")

all.censat.merge <- rbind(censat.merge,sst1.merge)


ggplot(all.censat.merge, aes(x=reorder(group, peaks), y=peaks))+geom_bar(stat="identity")+coord_flip()+ scale_y_break(c(50, 190))

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/All_SatRepeatsBar.pdf", 
       last_plot(), 
       height = 5,
       width=5)


