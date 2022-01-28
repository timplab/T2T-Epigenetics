
gcalls.gr <- readRDS("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls/whole_genome/HG002_nanonome_GpC_Smoothedbsojb.rds") 

coord=granges(gcalls.gr)

meth=data.frame(chr= seqnames(coord), start=start(coord), getMeth(gcalls.gr, type="smooth", what="perBase")) %>%
  mutate(score=rowMeans(.[,3:5])) %>%
  dplyr::select(chr,start,score) %>%
  filter(chr !="chrY") %>%
  mutate(end=start) %>%
  mutate(sites=1) %>%
  na.omit() %>%
  GRanges()

seqlevels(meth) <- seqlevels(BSgenome.t2t.v1.0.release)
seqinfo(meth) <- seqinfo(BSgenome.t2t.v1.0.release)
export.bw(meth, paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC.bw"))

freqmean=mean(meth$score)
freqsd=sd(meth$score)

zscore <- as.data.frame(meth) %>%
  mutate(z_score = (score - freqmean) / freqsd) %>%
  mutate(score=as.numeric(z_score)) %>%
  GRanges()

seqlevels(zscore) <- seqlevels(BSgenome.t2t.v1.0.release)
seqinfo(zscore) <- seqinfo(BSgenome.t2t.v1.0.release)
export.bw(zscore, paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_Smoothed_GC_zscore.bw"))

win=100000
chrs <- seqlevels(BSgenome.t2t.v1.0.release)

blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = chrs, width = win)

score1 <- coverage(meth, weight="score")
score2 <- coverage(meth, weight="sites")

binned_meth <- binnedSum(blocks, numvar = score1, "score") %>%
  as.data.frame() 
binned_sites <- binnedSum(blocks, numvar = score2, "sites") %>%
  as.data.frame() %>%
  dplyr::select(c(sites))


gc_bins <- cbind(binned_meth,binned_sites) %>%
  mutate(freq=score/sites) %>%
  na.omit()
freqmean=mean(gc_bins$freq)
freqsd=sd(gc_bins$freq)

gc_score <-  gc_bins %>%
  mutate(z_score = (freq - freqmean) / freqsd) %>%
  mutate(score=as.numeric(z_score)) %>%
  dplyr::select(-c(sites,freq,z_score)) %>%
  GRanges()

seqlevels(gc_score) <- seqlevels(BSgenome.t2t.v1.0.release)
seqinfo(gc_score) <- seqinfo(BSgenome.t2t.v1.0.release)

export.bw(gc_score, paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/mappability/HG002_100kbBinned_Smoothed_Accessibility.bw"))
