
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)
library(BSgenome.HG002.chrX)


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrM")

hg002_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_CpGmethylationFrequency_20kb.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrY") %>%
  filter(chromosome != "chrM")

censat.chm13 <- read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = c("chr", "start", "end", "name")) %>%
#  mutate(flag=case_when(chr != "chr9" | name != "HSAT3" ~ "keep", 
#                        TRUE ~ "remove")) %>%
#  filter(flag !="remove") %>%
  GRanges()

hg002_X <- read_tsv(paste0(dat, "/HG002/annotations/t2t_cenAnnotation.hg002_X.v1.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
  mutate(name = ifelse(grepl("Alu", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SATR", X4), "SATR", name)) %>%
  mutate(name = ifelse(grepl("ACRO", X4), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("GSATII", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("TAR", X4), "TAR", name)) %>%
  mutate(name = ifelse(grepl("TE", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MER", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MST", X4), "MST", name)) %>%
  mutate(name = ifelse(grepl("CER", X4), "CER", name)) %>%
  mutate(name = ifelse(grepl("L1", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SST", X4), "SST", name)) %>%
  mutate(name = ifelse(grepl("LSAU", X4), "LSAU", name)) %>%
  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("MSAT", X4), "MSAT", name)) %>%
  mutate(name = ifelse(grepl("novel", X4), "novel", name)) %>%
  mutate(name = ifelse(grepl("HERV", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("LT", X4), "TE", name)) %>%
  dplyr::select(c(X1, X2, X3, name)) %>%
  dplyr::rename("chr" =1, "start" = 2 ,"end" =3)

hg002_auto <- read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = c("chr", "start", "end", "name")) %>%
    filter(chr != "chrX")
  
censat.hg002 <- rbind(hg002_auto, hg002_X) %>%
 # mutate(flag=case_when(chr != "chr9" | name != "HSAT3" ~ "keep", 
 #                       TRUE ~ "remove")) %>%
#  filter(flag !="remove") %>%
    GRanges()

ggplot()+geom_histogram(data=as.data.frame(chm13_meth), aes(called_sites), fill="dodgerblue", alpha=.5)+geom_histogram(data=hg002_meth, aes(called_sites), fill="purple", alpha=.5)+scale_x_log10()


ggsave(
  paste0(figs, "CpG_Coverage_chm13vshg002.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 5
)



flags.hg002 <- hg002_meth %>% 
  filter(called_sites < 10 | called_sites > 100) %>%
  GRanges()

flags.chm13 <- chm13_meth %>% 
  filter(called_sites < 10 | called_sites > 100) %>%
  GRanges()

flags.hg002.censat <- FindOvls(flags.hg002, censat.hg002)
length(flags.hg002.censat$seqnames)/length(as.data.frame(flags.hg002)$seqnames)
flags.chm13.censat <- FindOvls(flags.chm13, censat.chm13)
length(flags.chm13.censat$seqnames)/length(as.data.frame(flags.chm13)$seqnames)

SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "gap-rDNA")


t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")


hg002.auto <- as.data.frame(t2t.cpg.loci[seqnames(t2t.cpg.loci) != "chrX"])

loci.hg002 <- findLoci(pattern = "CG",
                       subject = BSgenome.HG002.chrX::BSgenome.HG002.chrX, 
                       strand = "+")

loci.hg002 <- as.data.frame(loci.hg002)

t2t.loci.hg002 <- rbind(hg002.auto,loci.hg002) %>%
  GRanges()

total_cpgs_hg002 <- FindOvls(t2t.loci.hg002, censat.hg002) %>%
  group_by(name) %>%
  summarise(sum=sum(width)) %>%
  filter(name %in% SAT) %>%
  distinct() %>%
  rename(sum = "hg002_sum")

total_cpgs_chm13 <- FindOvls(t2t.cpg.loci, censat.chm13) %>%
  group_by(name, seqnames) %>%
  summarise(sum=sum(width)) %>%
  filter(name %in% SAT) %>%
  distinct()%>%
  rename(sum = "chm13_sum") 

flags.hg002.sum <- flags.hg002.censat %>%
  mutate(total=length(flags.hg002.censat$seqnames)) %>%
  group_by(seqnames, name) %>%
  summarise(bad_hg002_sum=sum(num_motifs_in_group), bad_hg002_frac=bad_hg002_sum/total) %>%
  filter(name %in% total_cpgs_hg002$name) %>%
  merge(total_cpgs_hg002, by="name") %>%
  mutate(frac_total = bad_hg002_sum/hg002_sum) %>%
  distinct()

flags.chm13.sum <- flags.chm13.censat %>%
  mutate(total=length(flags.chm13.censat$seqnames)) %>%
  group_by(seqnames, name) %>%
  summarise(bad_chm13_sum=sum(num_motifs_in_group), bad_chm13_frac=bad_chm13_sum/total) %>%
  filter(name %in% total_cpgs_chm13$name) %>%
  merge(total_cpgs_chm13, by="name") %>%
  mutate(frac_total = bad_chm13_sum/chm13_sum) %>%
  distinct()

flags.hg002 <- as.data.frame(hg002_meth) %>% 
  mutate(flags=case_when(called_sites < 10 | called_sites > 100 ~ "1", 
                         TRUE ~ "0")) %>%
  mutate(flags=as.numeric(flags)) %>%
  GRanges()

flags.chm13 <- as.data.frame(chm13_meth) %>% 
  mutate(flags=case_when(called_sites < 10 | called_sites > 100 ~ "1", 
                         TRUE ~ "0")) %>%
  mutate(flags=as.numeric(flags)) %>%
  GRanges()


blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = list, width = 10000)
score1 <- coverage(flags.chm13, weight="called_sites")
score2 <- coverage(flags.chm13, weight="flags")
binned_sites <- binnedSum(blocks, numvar = score1, "called_sites") %>%
  as.data.frame()
binned_flags <- binnedSum(blocks, numvar = score2, "flags") %>%
  as.data.frame()
chm13_flagged_bins <- merge(binned_sites, binned_flags, by=c("seqnames", "start","end")) %>%
  mutate(frac=flags/called_sites) %>%
  filter(frac > .9)

score1 <- coverage(flags.hg002, weight="called_sites")
score2 <- coverage(flags.hg002, weight="flags")
binned_sites <- binnedSum(blocks, numvar = score1, "called_sites") %>%
  as.data.frame()
binned_flags <- binnedSum(blocks, numvar = score2, "flags") %>%
  as.data.frame()
hg002_flagged_bins <- merge(binned_sites, binned_flags, by=c("seqnames", "start","end")) %>%
  mutate(frac=flags/called_sites) %>%
  filter(frac > .9)

write.table(hg002_flagged_bins, paste0(dat, "/chm13_final_beds/HG002_flagged_CpG_CoverageBins.bed"), col.names = T, row.names = F, quote=F, sep = "\t")

write.table(chm13_flagged_bins, paste0(dat, "/chm13_final_beds/CHM13_flagged_CpG_CoverageBins.bed"), col.names = T, row.names = F, quote=F, sep = "\t")
