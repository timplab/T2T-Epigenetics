#!/usr/bin/Rscript

library(tidyverse)
path="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome"
files=list.files(path, pattern= "_bismark.out")

bis=data.frame()
for (file in files){
  dat <- read_tsv(paste0(path, "/", file), col_names=c("chr", "start", "strand", "meth", "unmeth", "motif", "sequence")) %>%
    mutate(run=file)
  bis = rbind(dat, bis)
}

bis.merged=bis  %>%
  group_by(start) %>%
  summarise(start=start, strand=strand,meth=sum(meth), unmeth=sum(unmeth), motif=motif,sequence=sequence) %>%
  distinct()

out="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled"

write.table(bis.merged, paste0(out, "/HG002_nanonome1_winnowmapk15_merged_GpC_bismark.out"),
          col.names = F, row.names = F, sep="\t",quote = FALSE)
