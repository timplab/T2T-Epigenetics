
library(rtracklayer)
smooth_gc <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC.bw")
hg002_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
  GRanges()

acros <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO/ACRO_peaks.bed", col_names=c("chr", "start", "end")) %>%
  mutate(ID=row_number()) %>%
  mutate(len=end-start) %>%
  mutate(start=start+(len/2), end=start) %>%
  mutate(start=start-1000, end=end+1000) %>%
  GRanges()

ovls <- FindOvls(smooth_gc,acros) 

ovls.stat <- ovls %>%
  group_by(ID) %>%
  mutate(minstart=min(start), maxend=max(end), cen=(maxend-minstart)/2) %>%
  mutate(start=start-minstart,end=end-minstart) %>%
  arrange(start) %>%
  group_by(start) %>%
  summarise(meth=median(score))


ggplot(ovls.stat, aes(x=start, y=meth))+labs(x= "Position", y = "GpC Methylation")+coord_cartesian(xlim=c(500,1500))+geom_smooth(method="loess", span=.1, se=F)+geom_vline(xintercept = 1000, linetype="dashed")

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/ACRO/ACRO_Peak_position.pdf", 
       last_plot(), 
       height = 5,
       width=5)



sst1 <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/SST1_CEN-array-mon_NONCEN-array_withLabels_PeakIntersections.bed", col_names=c("chr", "start", "end", "rep_chr", "rep_start", "rep_end", "rep_strand", "group")) %>%
  mutate(ID=row_number()) %>%
  mutate(len=end-start) %>%
  mutate(start=start+(len/2), end=start) %>%
  mutate(start=start-1000, end=end+1000) %>%
  GRanges()

ovls <- FindOvls(smooth_gc,sst1) 

ovls.stat <- ovls %>%
  group_by(ID) %>%
  mutate(minstart=min(start), maxend=max(end), cen=(maxend-minstart)/2) %>%
  mutate(start=start-minstart,end=end-minstart) %>%
  arrange(start) %>%
  group_by(start) %>%
  summarise(meth=median(score))


ggplot(ovls.stat, aes(x=start, y=meth))+geom_vline(xintercept = 2000, linetype="dashed")+labs(x= "Position", y = "GpC Methylation")+coord_cartesian(xlim=c(500,1500))+geom_smooth(method="loess", span=.1, se=F)+geom_vline(xintercept = 1000, linetype="dashed")

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/SST1/SST1_Peak_position.pdf", 
       last_plot(), 
       height = 5,
       width=5)
