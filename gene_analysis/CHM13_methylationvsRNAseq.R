library(tidyverse)
library(RColorBrewer)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

protein_coding <- read_tsv(paste0(dat,"/gene_expression/dictionary.tsv")) %>%
  filter(gene_biotype == "protein_coding")

regs.all <- read_delim(paste0(dat, "/revision_analysis/promoter_clustering/RNAseq/CHM13.combined.v4.protein_coding.tss_RNAseq-cgi.bed"), col_names=c("chr","start","end","strand","name","gene_start","gene_end", "cov", "fpkm", "tpm", "cgi_chr","cgi_start","cgi_end"),delim="\t")  %>%
  filter(name %in% protein_coding$gene_id) %>%
  mutate(quantile=case_when(fpkm == 0 ~ "no_expression", 
                            fpkm > 0 & fpkm < 15 ~ "medium",
                            fpkm >=15 ~ "high"))

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
  group_by(ID,methylated_frequency,name, fpkm,region) %>% 
  summarise(gb_meth = sum(called_sites_methylated)/sum(called_sites)) %>%
  group_by(ID,name,fpkm,region) %>%
  summarize(frac_meth=mean(methylated_frequency))



cgi <- ovls %>%
  filter(region=="cgi")
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(cgi$fpkm, cgi$frac_meth)
title <- sprintf("N = %d r = %.3f", nrow(cgi), c)
ggplot(cgi, aes(cgi$frac_meth, log2(cgi$fpkm))) +
  geom_bin2d(bins=75) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Methylation") +
  ylab("RNA-seq FPKM") +
  theme_bw(base_size=20) +
  ggtitle(title)

ggsave(paste0(figs, "/CGIvsRNAseq.pdf"),
       last_plot(), 
       height=5,
       width=7)



gb <- ovls %>%
  filter(region=="gb")
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(gb$fpkm, gb$frac_meth)
title <- sprintf("N = %d r = %.3f", nrow(gb), c)
ggplot(gb, aes(gb$frac_meth, log2(gb$fpkm))) +
  geom_bin2d(bins=75) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Methylation") +
  ylab("RNA-seq FPKM") +
  theme_bw(base_size=20) +
  ggtitle(title)

ggsave(paste0(figs, "/GBvsRNAseq.pdf"),
       last_plot(), 
       height=5,
       width=7)