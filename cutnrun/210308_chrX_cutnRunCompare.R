
library(tidyverse)
library(rtracklayer)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


name="_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig"
path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/210204_alignments/022021_CUTRUN"

bw <- paste0(path,"/H3K9me3/CHM13_H3K9me3",name)
rep1 <- import(bw, format="bigwig")

bw <- paste0(path,"/IgG/CHM13_IgG",name)
igg1 <- import(bw, format="bigwig")

name="_cutnrun_losalt.F3852.over.IlluminaPCRfree_v1.0-assembly_51mers_single_mrg_meryl_sort.bigwig"
path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/210204_alignments/122920_CUTRUN"

bw <- paste0(path,"/H3K9me3/CHM13_H3K9me3",name)
rep2 <- import(bw, format="bigwig")

bw <- paste0(path,"/IgG/CHM13_IgG",name)
igg2 <- import(bw, format="bigwig")



reps <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0_plus38Y_repeatmasker.out.bed"), col_names = F)
SAT = c("ALR/Alpha", "GSATII", "GSATX","HSAT4")
reps.sat <- reps %>%
  filter(X4 %in% SAT) %>%
  dplyr::rename("chr"=1, "start"=2, "end"=3, "rep.name"=4)%>%
  select(chr, start, end, rep.name) %>%
  filter(chr == "chrX") %>%
  mutate(rep.name = ifelse(grepl("GSAT", rep.name), "GSAT", rep.name)) 

censat.gr <- GRanges(reps.sat)

keepi <- findOverlaps(rep1,censat.gr)
freq.matched <- rep1[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

rep1.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="h3k9me3") %>%
  mutate(rep="rep1")


keepi <- findOverlaps(igg1,censat.gr)
freq.matched <- igg1[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

igg1.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="igg") %>%
  mutate(rep="rep1")


# rep 2

keepi <- findOverlaps(rep2,censat.gr)
freq.matched <- rep2[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

rep2.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="h3k9me3") %>%
  mutate(rep="rep2")

keepi <- findOverlaps(igg2,censat.gr)
freq.matched <- igg2[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

igg2.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="igg") %>%
  mutate(rep="rep2")

all.h3k9 <- rbind(rep1.censat, rep2.censat, igg1.censat,igg2.censat)


treat.total1 <- sum(rep1$score)
control.total1 <- sum(igg1$score)
ratio1=control.total1/treat.total1

treat.total2 <- sum(rep2$score)
control.total2 <- sum(igg2$score)
ratio2=control.total2/treat.total2

stat.sum <- all.h3k9 %>%
  group_by(rep, treat,rep.name) %>%
  summarise(reads=sum(score)) %>%
  spread(treat, reads) %>%
  ungroup() %>%
  group_by(rep.name) %>%
  summarise(enrich=if_else(rep=="rep1",log2((h3k9me3/igg)*ratio1),log((h3k9me3/igg)*ratio2)))

ggplot(stat.sum, aes(x=rep.name, y=enrich, fill=rep.name))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = 1.5)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")

ggsave(
  paste0(figs, "/", "CHM13_chrXONLY_rmannotation.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)

name="_cutnrun-202021_losalt.F3852.over.IlluminaPCRfree-chm13-v1-autosome_hg002-chrX_51mers_mrg-meryl_sort.bigwig"
path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/210301_alignments/022021_CUTRUN"

bw <- paste0(path,"/H3K9me3/HG002_H3K9me3",name)
rep1 <- import(bw, format="bigwig")

bw <- paste0(path,"/IgG/HG002_IgG",name)
igg1 <- import(bw, format="bigwig")




reps <- read_tsv(paste0(dat, "/HG002/annotations/HG002_v0.9.chrX_repeatmasker.out.bed"), col_names = F)
SAT = c("ALR/Alpha", "GSATII", "GSATX","HSAT4")
reps.sat <- reps %>%
  filter(X4 %in% SAT) %>%
  dplyr::rename("chr"=1, "start"=2, "end"=3, "rep.name"=4)%>%
  select(chr, start, end, rep.name) %>%
  filter(chr == "chrX") %>%
  mutate(rep.name = ifelse(grepl("GSAT", rep.name), "GSAT", rep.name)) 

censat.gr <- GRanges(reps.sat)

keepi <- findOverlaps(rep1,censat.gr)
freq.matched <- rep1[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

rep1.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="h3k9me3") %>%
  mutate(rep="rep1")


keepi <- findOverlaps(igg1,censat.gr)
freq.matched <- igg1[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

igg1.censat <- as.data.frame(freq.matched) %>%
  mutate(treat="igg") %>%
  mutate(rep="rep1")


all.h3k9 <- rbind(rep1.censat, igg1.censat)


treat.total1 <- sum(rep1$score)
control.total1 <- sum(igg1$score)
ratio1=control.total1/treat.total1


stat.sum <- all.h3k9 %>%
  group_by(rep, treat,rep.name) %>%
  summarise(reads=sum(score)) %>%
  spread(treat, reads) %>%
  ungroup() %>%
  group_by(rep.name) %>%
  summarise(enrich=log2((h3k9me3/igg)*ratio1))

ggplot(stat.sum, aes(x=rep.name, y=enrich, fill=rep.name))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = 1.5)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")

ggsave(
  paste0(figs, "/", "HG002_chrXONLY_rmannotation.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)