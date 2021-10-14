library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(GenomicRanges)
library(zoo)
# load meth 

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated)



iso <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/transcripts.meth.and.cutandrun.bed") %>%
  rename( "chr"=`#chr`)

novel <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.Extra_paralogTrue.bed", col_names = c("chr", "start", "end", "gene"))

reps <- merge(iso, novel, by = c("chr", "start", "end", "gene")) %>%
  mutate(group=case_when(n_transctripts < 1 ~ "unexpressed",
                         TRUE ~ "expressed")) %>%
  dplyr::rename("direction"=strand )

flankn <- 3000
bodylen <- 1

# overlap with methylation calls and normalize by size 
l1_regions <- reps %>%
  mutate(start = start - flankn, end = end + flankn) %>%
  GRanges()

ovl <- findOverlaps(GRanges(chm13_meth), l1_regions)
genes.ovl <- as.data.frame(reps)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(GRanges(chm13_meth)[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  dplyr::mutate(dist = ifelse(direction == "+",start - gene_start, gene_end - start),
                dist = ifelse(dist < 0, dist/flankn,
                              ifelse(dist < genewidth,
                                     bodylen * dist / genewidth,
                                     bodylen + (dist - genewidth)/flankn)), 
                dist = round(dist,3)
  )

n_windows=500
chm13.ovl$cut = cut(chm13.ovl$dist, breaks=n_windows)

chm13.ovl.labs <- chm13.ovl %>%
  group_by(cut,group) %>%
  summarise(med = median(methylated_frequency), top = quantile(methylated_frequency, 0.75), bot = quantile(methylated_frequency, 0.25), n_genes = length(methylated_frequency)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  ungroup() %>%
  distinct() %>%
  arrange(min,group) %>%
  mutate(med_smooth = rollmean(med, 5, NA))


nums <- reps %>%
  mutate(total=n()) %>%
  group_by(group) %>%
  summarise(n=n(),total=total) %>%
  distinct()

# plot frequency with loess --- if doing SINEs don't use loess unless u want to die 
freqplot <- ggplot(chm13.ovl.labs,aes(x=min,y=med, color=group))+geom_line()+theme_classic()+
  annotate("text", x=1, y=1, label= paste0("expressed =", nums$n[1])) + 
  annotate("text", x=1, y=.9, label= paste0("unexpressed =", nums$n[2]))+theme(legend.position = "bottom", legend.direction="horizontal",axis.title.x=element_blank())+ylim(0,1)+scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1), labels = c(paste0("-",flankn/1e3,"kb"),"TTS","TES",paste0("+",flankn/1e3,"kb")))


box <- ggplot(reps,aes(x= "", y=log2(cutandrun_max), color=group))+geom_boxplot()+theme_classic()+theme(legend.position = "none")+labs(x="", y= "Log2(H3K4me2)")

plot_grid(freqplot,box,rel_widths = c(3, 1))


ggsave(
  paste0(dat, "/revision_analysis/gene_analysis/figures/NovelgeneExpression.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 12,
  height = 4,
)
