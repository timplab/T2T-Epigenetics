library(tidyverse)
library(zoo)

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
  mutate(cell_line = "HG03492")%>%
  mutate(norm_start = start - min(start))


meth_dat <- rbind(HG01109,HG01243,HG03098,HG03492)
sub <- meth_dat[seq(1, nrow(meth_dat), 100), ]
ggplot(data=sub, aes(x=norm_start/1e6,y=smoothed))+geom_line(stat="identity")+facet_wrap(~cell_line, ncol=1)+theme_classic()+labs(y="Methylation Frequency", x="[Mb]")
