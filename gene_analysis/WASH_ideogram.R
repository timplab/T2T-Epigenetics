library(BSgenome.t2t.v1)
library(tidyverse)
library(karyoploteR)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/sd/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


sample_cns <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/WASH.transcripts.and.meth.bed", skip =1, col_names = c("chr", "start", "end", "gene_id", "strand",  "n_transctripts",  "total_calls",  "meth_calls",  "frac_meth")) %>%
  mutate(ID = row_number()) %>% 
  mutate(n_transctripts, quantile_rank = ntile(n_transctripts,3)) %>%
  GRanges()

seqlevels(BSgenome.t2t.v1)
sample_cns <- sort(sample_cns)

pdf(paste0(figs, "/WASH_ideogram.pdf"), width = 5, height = 10)
kp <- plotKaryotype(BSgenome.t2t.v1,chromosomes=c("autosomal"))
kpAddBaseNumbers(kp)
#cn_col <- ifelse(sample_cns$cn>2, "red", "blue")
kpPlotMarkers(kp, data=sample_cns, labels=sample_cns$gene_id, text.orientation = "horizontal", r1=0.7, cex=0.8, marker.parts = c(0.2, 0.7, 0.1))
dev.off()