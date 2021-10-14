library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(rtracklayer)
library(bsseq)

t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/sd/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv"))
#chm13_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_CpGmethylationFrequency.tsv")) %>%
 # GRanges()

chm13.gr <- GRanges(chm13_meth)

h3k4me2 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/122920_CUTRUN/H3K4me2/CHM13_H3K4me2_cutnrun_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

h3k27me3 <- import("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/022021_CUTRUN/H3K27me3/CHM13_H3K27me3_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig")

genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/transcripts.meth.and.cutandrun.bed", skip =1, col_names = c("chr", "start", "end", "gene_id", "direction",  "n_transctripts",  "total_calls",  "meth_calls",  "frac_meth", "cutandrunmax")) %>%
  mutate(ID = row_number()) %>% 
  mutate(n_transctripts, quantile_rank = ntile(n_transctripts,3)) %>%
  GRanges()

meth.dat <- FindOvls(chm13.gr,genes)
h3k4m2.dat <- FindOvls(h3k4me2,genes)
h3k27me3.dat <- FindOvls(h3k27me3,genes)
cpg.dat <- FindOvls(t2t.cpg.loci,genes)

flankn <- 1e3
bodylen <- 3
SD_regions <- as.data.frame(genes) %>%
  mutate(start = start - flankn,
         end = end + flankn) %>%
  GRanges()

ovl <- findOverlaps(chm13.gr, SD_regions)
genes.ovl <- as.data.frame(genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(chm13.gr[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(direction == "-", gene_end - start, start - gene_start),
                                        dist = ifelse(dist < 0, dist/flankn,
                                                      ifelse(dist < genewidth,
                                                             bodylen * dist / genewidth,
                                                             bodylen + (dist - genewidth)/flankn)), 
                                        dist = round(dist,3)
  )


chm13.ovl.labs <- chm13.ovl %>%
  # mutate(reg = ifelse(Name %in% CT.genes$Name, "CT", "non-CT")) %>%
  group_by(quantile_rank,dist) %>%
  summarise(med_meth = median(methylated_frequency)) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = med_meth, color = as.factor(quantile_rank)))+
  geom_smooth(method = "loess", span = .015,se = F)+
 #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Aggregated Methylation Frequency") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()+ylim(0,1)#+facet_wrap(~gene_id)
p

ggsave(
  paste0(figs, "/allGenes.pdf"),
  plot = p,
  width = 8,
  height = 5
)


wash <- chm13.ovl %>%
  filter(grepl("NBPF", gene_id))

p <- ggplot(wash,aes( x = dist, y = methylated_frequency, color = as.factor(quantile_rank)))+
  geom_smooth(method = "loess", span = .3,se = F)+
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  #geom_point(alpha=.5,size = .5)+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()+ylim(0,1)#+facet_wrap(~n_transctripts+gene_id, scales = "free")
p

p <- ggplot(wash,aes( x = dist, y = methylated_frequency, color = as.factor(ID)))+
  geom_smooth(method = "loess", span = .3,se = F)+
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
    geom_point(alpha=.5,size = .5)+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()+ylim(0,1)+facet_wrap(~n_transctripts+gene_id, scales = "free")
p

cutnrun <- wash %>%
  dplyr::select(c(ID,quantile_rank,cutandrunmax)) %>%
  distinct()
box <- ggplot(cutnrun,aes(x= quantile_rank, y=log2(cutandrunmax), color=as.factor(quantile_rank)))+geom_boxplot()+theme_classic()+theme(legend.position = "none")+labs(x="", y= "Log2(H3K4me2)")

ggsave(
  paste0(figs, "/NBPF_Genes.pdf"),
  plot = p,
  width = 15,
  height = 15
)

##### H3K4me2 #############
ovl <- findOverlaps(h3k4me2, SD_regions)
genes.ovl <- as.data.frame(genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(h3k4me2[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(direction == "-", gene_end - start, start - gene_start),
         dist = ifelse(dist < 0, dist/flankn,
                       ifelse(dist < genewidth,
                              bodylen * dist / genewidth,
                              bodylen + (dist - genewidth)/flankn)), 
         dist = round(dist,3)
  )


chm13.ovl.labs <- chm13.ovl %>%
  # mutate(reg = ifelse(Name %in% CT.genes$Name, "CT", "non-CT")) %>%
  group_by(quantile_rank,dist) %>%
  summarise(med_meth = median(score)) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = log2(med_meth), color = as.factor(quantile_rank)))+
  geom_smooth(method = "loess", span = .015,se = F)+
  #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Log2(H3K4me2)") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()#+facet_wrap(~gene_id)
p

ggsave(
  paste0(figs, "/allGenesH3k4me2.pdf"),
  plot = p,
  width = 8,
  height = 5
)


##### H3k27me3 #############
ovl <- findOverlaps(h3k27me3, SD_regions)
genes.ovl <- as.data.frame(genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(h3k27me3[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(direction == "-", gene_end - start, start - gene_start),
         dist = ifelse(dist < 0, dist/flankn,
                       ifelse(dist < genewidth,
                              bodylen * dist / genewidth,
                              bodylen + (dist - genewidth)/flankn)), 
         dist = round(dist,3)
  )


chm13.ovl.labs <- chm13.ovl %>%
  # mutate(reg = ifelse(Name %in% CT.genes$Name, "CT", "non-CT")) %>%
  group_by(quantile_rank,dist) %>%
  summarise(med_meth = median(score)) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = log2(med_meth), color = as.factor(quantile_rank)))+
  geom_smooth(method = "loess", span = .02,se = F)+
  #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Log2(H3k27me3)") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()#+facet_wrap(~gene_id)
p

ggsave(
  paste0(figs, "/allGenesH3k27me3.pdf"),
  plot = p,
  width = 8,
  height = 5
)

###### CpG ##########
SD_regions <- as.data.frame(SD_regions) %>%
  group_by(quantile_rank) %>%
  mutate(total_bases = sum(width)) %>%
  GRanges()

ovl <- findOverlaps(GRanges(t2t.cpg.loci), SD_regions)
genes.ovl <- as.data.frame(genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end)

chm13.ovl <- as.data.frame(t2t.cpg.loci[queryHits(ovl),]) %>%
  bind_cols(genes.ovl) %>%
  dplyr::rename(seqnames = 1) %>%
  mutate(dist = ifelse(direction == "-", gene_end - start, start - gene_start),
         dist = ifelse(dist < 0, dist/flankn,
                       ifelse(dist < genewidth,
                              bodylen * dist / genewidth,
                              bodylen + (dist - genewidth)/flankn)), 
         dist = round(dist,3)
  )


total.bases = as.data.frame(SD_regions) %>%
  dplyr::select(c(quantile_rank,total_bases)) %>%
  distinct()

chm13.ovl.labs <- chm13.ovl %>%
  merge(total.bases, by = 'quantile_rank') %>%
  # mutate(reg = ifelse(Name %in% CT.genes$Name, "CT", "non-CT")) %>%
  group_by(quantile_rank,dist) %>%
  summarise(n = n()/total_bases) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = n, color = as.factor(quantile_rank)))+
  geom_smooth(method = "loess", span = .015,se = F)+
  #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Aggregated Methylation Frequency") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()#+facet_wrap(~gene_id)
p



p <- ggplot(chm13.ovl.labs,aes( x = dist, y = 1, fill = n))+
  geom_tile()+
  #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Aggregated Methylation Frequency") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()+facet_wrap(~quantile_rank, ncol=1)
p

ggsave(
  paste0(figs, "/allGenesCpG.pdf"),
  plot = p,
  width = 8,
  height = 5
)