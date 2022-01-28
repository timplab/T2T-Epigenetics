library(tidyverse)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg38)
t2t.cpg.loci <- findLoci(pattern = "CG",
                         subject = BSgenome.t2t.v1.0.release::BSgenome.t2t.v1.0.release, 
                         strand = "+")
length(t2t.cpg.loci)

g002.nanopore <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled_50kb.tsv") %>%
  mutate(method="nanopore")

hg002.bsseq <- read_tsv("/uru/Data/old_atium/Data/Nanopore/projects/ONT_HG002_bisulfite/bsseq/bismark/CpGCalls/004_0111_001_R1_val_1_bismark_bt2_pe.bismark.cov.gz", col_names = c("chromosome", "start", "end", "methylated_frequency", "called_sites_methylated", "called_sites_unmethylated")) %>%
  mutate(start=start-1, end=end-1,methylated_frequency=methylated_frequency/100, called_sites =called_sites_methylated+called_sites_unmethylated, num_motifs_in_group=1, group_sequence="CG",method="bsseq") %>%
  dplyr::select(chromosome, start, end,num_motifs_in_group,called_sites,called_sites_methylated, methylated_frequency,group_sequence, method) %>%
  GRanges()

chm13_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/methylation_frequency_50kb_split.tsv") %>%
  mutate(method="chm13_nanopore")

bsseq.ref <- FindOvls(t2t.cpg.loci,hg002.bsseq) %>%
  select(-c(width,strand)) %>%
  dplyr::rename("chromosome"=seqnames)

total.bsseq=sum(bsseq.ref$called_sites)
q1.bsseq=quantile(bsseq.ref$called_sites, .05) 
q2.bsseq=quantile(bsseq.ref$called_sites, .95) 

q1.nanopore=quantile(g002.nanopore$called_sites, .05) 
q2.nanopore=quantile(g002.nanopore$called_sites, .95) 
total.nanopore=sum(g002.nanopore$called_sites)

q1.chm13nanopore=quantile(chm13_meth$called_sites, .05) 
q2.chm13nanopore=quantile(chm13_meth$called_sites, .95) 
total.chm13nanopore=sum(chm13_meth$called_sites)

hg002.nanopore <- g002.nanopore %>%
  mutate(qual=case_when(called_sites < q1.nanopore ~ "bad", 
                        called_sites > q2.nanopore ~ "bad", 
                        TRUE ~ "good"))

chm13.nanopore <- chm13_meth %>%
  mutate(qual=case_when(called_sites < q1.chm13nanopore ~ "bad", 
                        called_sites > q2.chm13nanopore ~ "bad", 
                        TRUE ~ "good"))

bsseq.ref <- bsseq.ref %>%
  mutate(qual=case_when(called_sites < q1.bsseq ~ "bad", 
                        called_sites > q2.bsseq ~ "bad", 
                        TRUE ~ "good"))

rep.widths  <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_RepeatMaskerV2.bed") %>%
  mutate(width=chromEnd-chromStart) %>%
  filter(repClass %in% list2) %>%
 # filter(name != "ALR/Alpha") %>%
 # filter(name != "HSATI") %>%
  dplyr::rename("chrom"=`#chrom`, "start"=chromStart,"end"=chromEnd) %>%
  GRanges()
  

list=c("HOR", 'DHOR', "HSAT1", "HSAT2", "HSAT3", "HSAT4", "MON", "SST1_Composite", "ACRO", "GSAT")

censat.widths <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
mutate(name = ifelse(grepl("hsat1", name), "HSAT1", name)) %>%
mutate(name = ifelse(grepl("hsat2", name), "HSAT2", name)) %>%
mutate(name = ifelse(grepl("hsat3", name), "HSAT3", name)) %>%
mutate(name = ifelse(grepl("hsat4", name), "HSAT4", name)) %>%
mutate(name = ifelse(grepl("hsat5", name), "HSAT5", name)) %>%
mutate(name = ifelse(grepl("hsat6", name), "HSAT6", name)) %>%
mutate(name = ifelse(grepl("ct", name), "CT", name)) %>%
mutate(name = ifelse(grepl("bsat", name), "BSAT", name)) %>%
mutate(name = ifelse(grepl("dhor", name), "DHOR", name)) %>%
mutate(name = ifelse(grepl("hor", name), "HOR", name)) %>%
mutate(name = ifelse(grepl("mon", name), "MON", name)) %>%
mutate(name = ifelse(grepl("GSAT", name), "GSAT", name)) %>%
mutate(name = ifelse(grepl("ACRO", name), "ACRO", name)) %>%
mutate(name = ifelse(grepl("SST1_Composite", name), "SST1_Composite", name)) %>%
mutate(width=chromEnd-chromStart) %>%
filter(name %in% list) %>%
dplyr::rename("chrom"=`#chrom`, "start"=chromStart,"end"=chromEnd) %>%
GRanges()

strand(t2t.cpg.loci)="*"

cg.reps <- FindOvls(t2t.cpg.loci, censat.widths) %>%
  group_by(name) %>%
  summarise(total_cgs=n())

nanopore.reps <- FindOvls(GRanges(hg002.nanopore), censat.widths) %>%
  filter(qual=="good") %>%
  group_by(name) %>%
  summarise(n=n()) %>%
  mutate(method="nanopore")

chm13.nanopore.reps <- FindOvls(GRanges(chm13.nanopore), censat.widths) %>%
  filter(qual=="good") %>%
  group_by(name) %>%
  summarise(n=n()) %>%
  mutate(method="chm13_nanopore") 

bsseq.reps <- FindOvls(GRanges(bsseq.ref), censat.widths)%>%
  filter(qual=="good") %>%
  group_by(name) %>%
  summarise(n=n()) %>%
  mutate(method="bsseq") 

rep.all <- rbind(bsseq.reps,nanopore.reps,chm13.nanopore.reps) %>%
  merge(cg.reps) %>%
  group_by(name, method) %>%
  mutate(frac=n/total_cgs)


ggplot(rep.all, aes(x=name, y=frac, fill=method))+geom_bar(stat = "identity", position="dodge")

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"

ggsave(
  paste0(figs, "/FractionHighQuality_CENSAT.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)



strand(t2t.cpg.loci)="*"

cg.reps <- FindOvls(t2t.cpg.loci, censat.widths) %>%
  filter(name=="HOR") %>%
  group_by(name,seqnames) %>%
  summarise(total_cgs=n())

nanopore.reps <- FindOvls(GRanges(hg002.nanopore), censat.widths) %>%
  filter(qual=="good") %>%
  filter(name=="HOR") %>%
  group_by(name,seqnames) %>%
  summarise(n=n()) %>%
  mutate(method="nanopore")

chm13.nanopore.reps <- FindOvls(GRanges(chm13.nanopore), censat.widths) %>%
  filter(qual=="good") %>%
  filter(name=="HOR") %>%
group_by(name,seqnames) %>%
  summarise(n=n()) %>%
  mutate(method="chm13_nanopore") 

bsseq.reps <- FindOvls(GRanges(bsseq.ref), censat.widths)%>%
  filter(qual=="good") %>%
  filter(name=="HOR") %>%
group_by(name,seqnames) %>%
  summarise(n=n()) %>%
  mutate(method="bsseq") 

rep.all <- rbind(bsseq.reps,nanopore.reps,chm13.nanopore.reps) %>%
  merge(cg.reps) %>%
  group_by(name, method,seqnames) %>%
  mutate(frac=n/total_cgs) %>%
  filter(method != "bsseq")


ggplot(rep.all, aes(x=seqnames, y=frac, fill=method))+geom_bar(stat = "identity", position="dodge")

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/figures"

ggsave(
  paste0(figs, "/FractionHighQuality_HORonly.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 10,
  height = 5,
)

