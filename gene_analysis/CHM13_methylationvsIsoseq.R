library(tidyverse)
library(RColorBrewer)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

protein_coding <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/mitchell_merged_bed/protein_coding_names.txt"),col_names="gene")

regs.all <- read_delim(paste0(dat, "/revision_analysis/promoter_clustering/mitchell_merged_bed/chm13_v1_tss.CGI_ntranscripts_intersections.bed"), col_names=c("chr","start","end","strand","name","gene_start","gene_end", "transcripts", "cgi_chr","cgi_start","cgi_end"),delim="\t")  %>%
  filter(name %in% protein_coding$gene) %>%
  mutate(quantile=case_when(transcripts == 0 ~ "no_expression", 
                            transcripts > 0 & transcripts < 113 ~ "medium",
                            transcripts >=100 ~ "high")) 

regs.cgi <- regs.all  %>%
  mutate(region="cgi") 

regs.gb <- regs.all  %>%
  mutate(start=ifelse(strand=="+", cgi_end,gene_end))%>%
  mutate(end=ifelse(strand=="+", gene_end,cgi_end)) %>%
  mutate(len=end-start) %>%
  filter(len > 100) %>%
  select(-c(len))  %>%
  mutate(region="gb")

regs <- rbind(regs.gb,regs.cgi) %>%
  distinct() %>%
  mutate(ID=row_number()) %>%
  GRanges()

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated)

ovls <- FindOvls(GRanges(chm13_meth),regs) %>%
  group_by(ID,methylated_frequency,name, transcripts,region) %>% 
  summarise(gb_meth = sum(called_sites_methylated)/sum(called_sites)) %>%
  group_by(ID,name,transcripts,region) %>%
  summarize(frac_meth=mean(methylated_frequency))



cgi <- ovls %>%
  filter(region=="cgi")
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(cgi$transcripts, cgi$frac_meth)
title <- sprintf("N = %d r = %.3f", nrow(cgi), c)
ggplot(cgi, aes(cgi$frac_meth, log2(cgi$transcripts))) +
  geom_bin2d(bins=75) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Methylation") +
  ylab("Iso-seq") +
  theme_bw(base_size=20) +
  ggtitle(title)



ggsave(paste0(figs, "/GBvsIsoseq.pdf"),
       last_plot(), 
       height=5,
       width=7)


gb <- ovls %>%
  filter(region=="gb")
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(gb$transcripts, gb$frac_meth)
title <- sprintf("N = %d r = %.3f", nrow(gb), c)
ggplot(gb, aes(gb$frac_meth, log2(gb$transcripts))) +
  geom_bin2d(bins=75) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Methylation") +
  ylab("Iso-seq") +
  theme_bw(base_size=20) +
  ggtitle(title)


ggsave(paste0(figs, "/CGIvsIsoseq.pdf"),
       last_plot(), 
       height=5,
       width=7)