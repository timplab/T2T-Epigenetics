library('liftOver')
library('tidyverse')

path="/uru/Data/Nanopore/projects/nanonome/backup/haplotypes/"
output.dir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/GM12878_liftover"
gm.hg38.gr <- read_tsv(paste0(path,"GM12878_nanoNOMe.hap1.gpc.mfreq.txt.gz"), col_names = c("chr", "start", "strand", "Methylated", "Unmethylated")) %>%
  mutate(end=start) %>%
  GRanges()

t2t.ch <- import.chain("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/bisulfite/liftover/hg38.t2t-chm13-v1.0.over.chain")

chm.hg38 <- liftOver(gm.hg38.gr, t2t.ch)
chm.hg38.gr <- unlist(chm.hg38)
saveRDS(chm.hg38.gr,file=file.path(output.dir,"gm12878.gpcmethhap1.chm13.gr.rds"))
