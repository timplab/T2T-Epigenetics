library(tidyverse)
library(rtracklayer)
library(BSgenome.t2t.v1.0.release)
library(Repitools)
library(GenomicRanges)

meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/methylation_frequency_50kb_split.tsv") %>%
  group_by(chromosome) %>%
  mutate(end=ifelse(start==end, start+1, end)) %>%
  filter(chromosome != "chrY") %>%
  filter(chromosome != "chrM") %>%
  dplyr::select(chromosome, start,end,methylated_frequency) %>%
  dplyr::rename("score"=methylated_frequency) %>%
  GRanges()


seqlevels(meth) <- seqlevels(BSgenome.t2t.v1.0.release)[1:23]
seqdat <- as.data.frame(seqinfo(BSgenome.t2t.v1.0.release))[1:23,]

seqinfo(meth) <- Seqinfo(seqnames=rownames(seqdat),
                         seqlengths=seqdat$seqlengths,
                         isCircular=seqdat$isCircular,
                         genome=seqdat$genome)

export.bw(meth, paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/bigwigs/CHM13_CpG_methylationFrequency_50kb.bigwig"))
