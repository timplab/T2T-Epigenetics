library(BSgenome.t2t.v1)
library(tidyverse)
library(karyoploteR)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


synteny <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/GRCh38_Synteny_25Kb.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
  GRanges()

censat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
  dplyr::rename(chrom=`#chrom`, start=chromStart, end=chromEnd) %>%
  GRanges()

sample_cns <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/NBPF.transcripts.meth.and.cutandrunMerged.bed",col_names = c("chr", "start", "end", "gene_id", "score")) %>%
  mutate(ID = row_number()) %>% 
  GRanges()
chrom=c("chr1")
seqlevels(BSgenome.t2t.v1)
sample_cns <- sort(sample_cns)

pdf(paste0(figs, "/NBPF_ideogram.pdf"), width = 5, height = 10)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes=chrom,plot.type=2)
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kpPlotRegions(kp, data=sample_cns,data.panel=1,r0=.1, r1=1)
kpRect(kp, data= censat, y0=0, y1=1, col= colors, data.panel="ideogram", border=NA)
kpRect(kp, synteny,data.panel=2, y0=.4, y1=.6, r0=.4, r1=.6,col="orange")
dev.off()

#chr1:144,643,857-148,942,106
#detail.region <- toGRanges(data.frame("chr1", 144643857, 148942106))
pdf(paste0(figs, "/NBPF_ideogram.pdf"), width = 5, height = 10)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes=chrom,plot.type=2) #,zoom=detail.region)
kpAddBaseNumbers(kp)
kpPlotMarkers(kp, data=sample_cns,data.panel=1,r0=.1, r1=1,labels=sample_cns$gene_id,text.orientation = "horizontal",label.dist = 0.01, max.iter = 1000)
dev.off()
