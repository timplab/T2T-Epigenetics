library(tidyverse)

des <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDeserts.bed") %>%
  mutate(len=chromEnd-chromStart) %>%
  group_by(len) %>%
  summarise(sum = sum(len)) %>%
  mutate(desfrac=sum/217619915,genfrac=sum/3056916522) %>%
  arrange(len) %>%
  mutate(totdes=cumsum(desfrac),totgen=cumsum(genfrac)) 

perc <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_MarkerDeserts.bed") %>%
  mutate(len=chromEnd-chromStart) %>%
  filter(len > 50000) %>%
  summarise(sum = sum(len)) %>%
  mutate(desfrac=sum/217619915,genfrac=(sum/3056916522)*100)

ggplot(des, aes(x=len/1e3,y=totgen*100))+geom_line()+geom_point()+geom_vline(xintercept = 50000/1e3, linetype="dashed")+scale_x_log10()+labs(y="Percentage of the Genome", x = "Desert Length (Kb)")+
  theme(text = element_text(size=18))+theme_classic() + annotate("text", x = 500000/1e3, y = 5, label = paste0(round(perc$genfrac,2),"%"))

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures/MarkerDesertGenomePercent.pdf",
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 5,
)
