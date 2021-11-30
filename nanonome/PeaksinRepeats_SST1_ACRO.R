library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(ggbreak)
options(scipen = 100)


sst1.peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/SST1_CEN-array-mon_NONCEN-array_withLabels_PeakIntersections.bed", col_names=c("chr", "peak_start", "peak_end", "rep_chr", "rep_start", "rep_end", "rep_strand", "group")) %>%
  group_by(group,chr) %>%
  summarize(n=n()) 

sst1.all <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/revision/SST1_CEN-array-mon_NONCEN-array_withLabels.txt")%>%
  mutate(width=end-start) %>%
  group_by(group,`#chrom`) %>%
  dplyr::rename("chr"=`#chrom`) %>%
  summarise(width=sum(width))

sst1.merge  <- merge(sst1.peaks,sst1.all, all=TRUE) %>%
 # filter(group != "mon_CEN") %>%
  mutate(peaks=(n/width)*1e4)

sst1.merge[is.na(sst1.merge)] <- 0

ggplot(sst1.merge, aes(x=reorder(chr, peaks), y=peaks, fill=group))+geom_bar(stat="identity", position="dodge")+coord_flip()+geom_text(stat='identity', aes(label=peaks), vjust=-1)

write.table(sst1.merge, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/SST1_table.txt", col.names = T, row.names = F, quote = F)

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/All_SatRepeatsBar.pdf", 
       last_plot(), 
       height = 5,
       width=5)



acro.peaks <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/ACRO_peaks_PeakIntersections.bed", col_names=c("chr", "peak_start", "peak_end", "rep_chr", "rep_start", "rep_end", "rep_type", "where", "width", "array_status")) %>%
  group_by(chr,array_status) %>%
  summarize(n=n()) 

acro.all <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO/ACRO.coords.bed", col_names=c("chr", "start", "end", "type", "where", "width", "array_status")) %>%
  group_by(chr,array_status) %>%
  summarize(width=sum(width))

acro.merge  <- merge(acro.peaks,acro.all, all=TRUE) %>%
  # filter(group != "mon_CEN") %>%
  mutate(peaks=(n/width)*1e5)

acro.merge[is.na(acro.merge)] <- 0

ggplot(acro.merge, aes(x=reorder(chr, peaks), y=peaks, fill=array_status))+geom_bar(stat="identity", position="dodge")+coord_flip()+geom_text(stat='identity', aes(label=peaks), vjust=-1)+facet_wrap(~array_status)

write.table(acro.merge, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO_table.txt", col.names = T, row.names = F, quote = F)

