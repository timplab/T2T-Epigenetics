library('tidyverse')
library('ggpubr')


genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/satellites/TE_locations/Alus_in_genes.bed", col_names = c("chr", "start", "end", "Alu", "Direction", "strand", "rep_start", "rep_end", "score2", "dir", "len", "num","HOT","gene_chr", 'gene_start',"gene_end", "name", "strand", "n_transcripts", "calls", "calls2", "frac_meth", "cutnrun_max")) %>%
  mutate(ID=row_number()) %>%
  GRanges()

dict <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/dictionary.tsv")
HG002.exp <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/stringtie2/HG002_SRR13086640/HG002_SRR13086640_abun.txt") %>%
  dplyr::rename("gene_id" = `Gene ID`, "gene_name"= `Gene Name`) %>%
  merge(dict, by="gene_id") %>%
  dplyr::rename("name"=gene_name.y) %>%
  dplyr::select(c(Coverage,FPKM,name))

HG002.genes <- merge(as.data.frame(genes),HG002.exp,by="name")%>%
  dplyr::select(c(seqnames, start,end,Alu,HOT,Coverage,FPKM,ID))%>%
  GRanges()


dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

hg002_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
  GRanges()

hg002_gc <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_GpCmethylationFrequency.tsv")) %>%
  GRanges()

meth.ovls <- FindOvls(HG002.genes,hg002_meth) %>%
  group_by(ID) %>%
  mutate(rep_meth=mean(methylated_frequency)) %>%
  mutate(expressed=case_when(FPKM > 1 ~ "on", 
                             TRUE ~ "off"))

ggplot(meth.ovls,aes(x=as.factor(expressed), y=rep_meth))+geom_violin()+geom_boxplot(outlier.shape = NA, width=.1)+stat_compare_means(method = "wilcox.test")


meth.ovls.stat <- meth.ovls %>%
  dplyr::select(expressed, rep_meth) %>%
  group_by(expressed) %>%
  summarise(med_meth=median(rep_meth))

gc.ovls <- FindOvls(HG002.genes,hg002_gc) %>%
  group_by(ID) %>%
  mutate(rep_meth=mean(methylated_frequency)) %>%
  mutate(expressed=case_when(FPKM > 1 ~ "on", 
                             TRUE ~ "off"))

ggplot(gc.ovls,aes(x=as.factor(expressed), y=rep_meth))+geom_violin()+geom_boxplot(outlier.shape = NA, width=.1)+stat_compare_means(method = "wilcox.test")

gc.ovls.stat <- gc.ovls %>%
  dplyr::select(expressed, rep_meth) %>%
  group_by(expressed) %>%
  summarise(med_meth=median(rep_meth))

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()
meth.ovls.chm13 <- FindOvls(genes,chm13_meth) %>%
  group_by(ID) %>%
  mutate(rep_meth=mean(methylated_frequency)) %>%
  mutate(expressed=case_when(n_transcripts > 3 ~ "on", 
                             TRUE ~ "off"))

ggplot(meth.ovls.chm13,aes(x=as.factor(expressed), y=rep_meth))+geom_violin()+geom_boxplot(outlier.shape = NA, width=.1)+stat_compare_means(method = "wilcox.test")
ggplot(meth.ovls.chm13,aes(expressed, fill=HOT))+geom_bar()



flankn <- 100
bodylen <- 500

# overlap with m
# overlap with methylation calls and normalize by size 
l1_regions <- HG002.genes %>%
  dplyr::mutate(gene_start = start, gene_end = end) %>%
  mutate(reg_end = ifelse(Direction == "+", end + flankn, start-flankn)) %>%
  mutate(reg_start = ifelse(Direction=="+", end - bodylen, start+bodylen)) %>%
  mutate(start=ifelse(Direction=="+", reg_start,reg_end)) %>%
  mutate(end=ifelse(Direction=="+", reg_end,reg_start)) %>%
  GRanges()

ovl <- findOverlaps(hg002_meth, l1_regions)
genes.ovl <- as.data.frame(HG002.genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(hg002_meth[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(Direction == "+",gene_end - start,start - gene_start),
                #dist = ifelse(dist < 0, dist/flankn,
                #              ifelse(dist < genewidth,
                #                     bodylen * dist / genewidth,
                #                     bodylen + (dist - genewidth)/flankn)), 
                dist = plyr::round_any(dist,5)*-1)
#  )

n_windows=1000
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group,escape) %>%
  summarise(med = median(methylated_frequency),ID=ID,escape=escape) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  group_by(group,escape) %>%
  arrange(min) %>%
  mutate(med_smooth = rollmean(med, 2, NA))

# plot frequency with loess --- if doing SINEs don't use loess unless u want to die 
freqplot <- ggplot(chm13.ovl.labs,aes(x=min,y=med_smooth, color=group))+theme(legend.position = "left", legend.direction="vertical",axis.title.x=element_blank())+ylim(0,1)+ geom_smooth(method="loess",se=F, span=.3)+facet_wrap(~escape)#+geom_line()

# plot heatmap
plot <- ggplot(chm13.ovl,aes(x=dist,y=-ID,fill=methylated_frequency))+ geom_tile()+scale_fill_gradient(low = "blue", high = "red", na.value = NA, limits=c(0,1))+theme(legend.position = "left", legend.direction="vertical")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())+facet_wrap(~group+escape,ncol=1)#+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))

# stack density and heatmap
plot_grid(freqplot, plot, ncol=1, align = "v", rel_heights=c(1/3,1))


#ggsave(
#  paste0(figs, "/", rep, "_Heatmap.pdf"),
#  plot = last_plot(),
#  scale = 1,
#  width = 8,
#  height = 12,
#)
#
#ggsave(
#  paste0(figs, "/" , rep, "_Heatmap.png"),
#  plot = last_plot(),
#  scale = 1,
#  width = 8,
#  height = 12,
#)
