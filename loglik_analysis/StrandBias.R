knitr::opts_chunk$set(echo = FALSE, warning=FALSE)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


censat.meth <- read_tsv(paste0(dat, "/methylation_calls/repeat_intersect/50kb_methylationCalls_censatV2Intersect.bed"), col_names = F)

censat.meth.stat <- censat.meth %>%
  filter(X6 < -1.5 | X6 > 1.5) %>%
  mutate(mcall=case_when(X6 < -1.5 ~ "0", 
                         X6 > 1.5 ~ "1")) %>%
  mutate(mcall=as.numeric(mcall)) %>%
  group_by(X1, X2,X4,X10) %>%
  summarize(cov=n(),methylated_frequency=mean(mcall))

SAT=c("HOR", "DHOR", "MON", "CT", "HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5")

censat.meth.sub <- censat.meth.stat %>%
  filter(X10 %in% SAT)
ggplot(censat.meth.sub, aes(x=X10, y=methylated_frequency, fill=X4))+geom_boxplot(outlier.shape = NA)

ggsave(
  paste0(figs, "/CensatMethylationStrand.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 6,
)

ggplot(censat.meth.sub, aes(x=X10, y=cov, fill=X4))+geom_boxplot(outlier.shape = NA)+ylim(0,45)

ggsave(
  paste0(figs, "/CensatCovStrand.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 6,
)


censat.meth.strand <- censat.meth.sub %>%
  select(-c(cov)) %>%
  group_by(X2) %>%
  spread(X4,methylated_frequency) %>%
  rename(X1="chr", X2="start", `-`="minus", `+`="plus") %>%
  filter(chr != "chr9" & X10 != "HSAT3")

pal <- wes_palette("Zissou1", 100, type = "continuous")
sub <- censat.meth.strand %>%
  filter(X10=="HOR")
ggplot(censat.meth.strand, aes(x=minus, y=plus))+geom_density_2d_filled(contour_var = "ndensity")+facet_wrap(~X10)


ggplot(sub, aes(x=minus, y=plus))+stat_bin_2d(bins=50)+facet_wrap(~X10)+scale_fill_gradientn(colours = pal) 

stat.censat <- censat.meth %>%
  filter(X6 < -1.5 | X6 > 1.5) %>%
  mutate(mcall=case_when(X6 < -1.5 ~ "0", 
                         X6 > 1.5 ~ "1")) %>%
  mutate(mcall=as.numeric(mcall)) %>%
  group_by(X1, X2,X10) %>%
  summarize(cov=n(),methylated_frequency=mean(mcall)) %>%
  group_by(X10) %>%
  summarize(coverage=mean(cov), methylation=mean(methylated_frequency))

write.table(stat.censat, file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/repeat_fasta/Average_coverage.csv",  quote = F, sep = ",", row.names = F,
            col.names = TRUE)
  