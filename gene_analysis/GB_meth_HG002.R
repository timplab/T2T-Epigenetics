library('tidyverse')
library('ggpubr')
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")


dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

hg002_meth <- read_tsv(paste0(dat, "/revision_analysis/HG002_pooled/HG002_CpG_methylationFrequency_pooled.tsv")) %>%
  GRanges()

hg002_gc <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/chm13_whole_genome/pooled/HG002_nanonome_GpCmethylationFrequency.tsv")) %>%
  GRanges()


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/figures"
protein_coding <- read_tsv(paste0(dat,"/revision_analysis/promoter_clustering/mitchell_merged_bed/protein_coding_names.txt"),col_names="gene")

genes <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/transcripts.meth.and.cutandrun.bed", skip =1, col_names = c("chr", "start", "end", "name", "direction",  "n_transctripts",  "total_calls",  "meth_calls",  "frac_meth", "cutandrunmax")) %>%
  mutate(ID = row_number())


dict <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/dictionary.tsv")
HG002.exp <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/gene_expression/stringtie2/HG002_SRR13086640/HG002_SRR13086640_abun.txt") %>%
  dplyr::rename("gene_id" = `Gene ID`, "gene_name"= `Gene Name`) %>%
  merge(dict, by="gene_id") %>%
  dplyr::rename("name"=gene_name.y) %>%
  dplyr::select(c(Coverage,FPKM,name,gene_biotype)) %>%
  filter(gene_biotype=="protein_coding")

HG002.genes <- merge(as.data.frame(genes),HG002.exp,by="name") %>%
  mutate(expressed=ntile(FPKM,5)) %>%
  GRanges()

flankn <- 1e3
bodylen <- 3
SD_regions <- as.data.frame(HG002.genes) %>%
  mutate(start = start - flankn,
         end = end + flankn) %>%
  GRanges()

ovl <- findOverlaps(hg002_meth, SD_regions)
genes.ovl <- as.data.frame(HG002.genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(hg002_meth[queryHits(ovl),]) %>%
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
  group_by(expressed,dist) %>%
  summarise(med_meth = median(methylated_frequency)) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = med_meth, color = as.factor(expressed)))+
  geom_smooth(method = "loess", span = .01,se = F)+
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
  paste0(figs, "/HG002_allGenes_CpG.pdf"),
  plot = p,
  width = 8,
  height = 5
)

ovl <- findOverlaps(hg002_gc, SD_regions)
genes.ovl <- as.data.frame(HG002.genes)[subjectHits(ovl),] %>%
  dplyr::mutate(genewidth = end - start) %>%
  dplyr::rename(gene_start = start, gene_end = end) 

chm13.ovl <- as.data.frame(hg002_gc[queryHits(ovl),]) %>%
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
  group_by(expressed,dist) %>%
  summarise(med_meth = median(methylated_frequency)) %>%
  distinct()

p <- ggplot(chm13.ovl.labs,aes( x = dist, y = med_meth, color = as.factor(expressed)))+
  geom_smooth(method = "loess", span = .01,se = F)+
  #  geom_point()+
  #  geom_smooth(se=T)+
  # geom_vline(xintercept = 0) +
  # geom_vline(xintercept = bodylen) +
  scale_x_continuous(breaks= c(-1,0,bodylen,bodylen + 1),
                     labels = c(paste0("-",flankn/1e3,"kb"),"Start","End",paste0("+",flankn/1e3,"kb"))) +
  labs( x = "Genomic Position", y = "Aggregated Methylation Frequency") +
  theme(legend.background = element_rect(color = "black")) +
  theme_classic()+ylim(.1,.75)#+facet_wrap(~gene_id)
p
ggsave(
  paste0(figs, "/HG002_allGenes_GpC.pdf"),
  plot = p,
  width = 8,
  height = 5
)

tss <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/gene_analysis/CHM13.combined.v4.TSSonly.bed", col_names=c("chr", "start", "end","strand", "gene_name", "gene_biotype"), delim="  ") %>%
  separate(gene_name, c("gene_name", "transcript"), "-") %>%
  group_by(gene_name) %>%
  filter(max(end)-min(start) < 1000) %>%
  select(c(gene_name)) %>%
  distinct()



protein_coding <- read_tsv(paste0(dat,"/gene_expression/dictionary.tsv")) %>%
  filter(gene_biotype == "protein_coding") %>%
  dplyr::rename("name"=gene_name)
#filter(gene_name %in% tss$gene_name)

regs.all <- read_delim(paste0(dat, "/revision_analysis/promoter_clustering/RNAseq/HG002/HG002.combined.v4.protein_coding.tss_RNAseq-cgi.bed"), col_names=c("chr","start","end","strand","gene_id","gene_start","gene_end", "cov", "fpkm", "tpm", "cgi_chr","cgi_start","cgi_end"),delim="\t")  %>%
  merge(protein_coding, by="gene_id") %>%
  # mutate(quantile=ntile(fpkm,4)) %>%
  mutate(quantile=ntile(fpkm,5))


regs.cgi <- regs.all  %>%
  mutate(region="cgi") 

regs.gb <- regs.all  %>%
  mutate(start=ifelse(strand=="+", cgi_end+2000,gene_end))%>%
  mutate(end=ifelse(strand=="+", gene_end,cgi_end-2000)) %>%
  mutate(len=end-start) %>%
  filter(len > 1000) %>%
  dplyr::select(-c(len))  %>%
  mutate(region="gb")


regs <- rbind(regs.gb,regs.cgi) %>%
  distinct() %>%
  # mutate(start=ifelse(strand=="+", gene_start,gene_end))%>%
  #  mutate(end=ifelse(strand=="+", gene_end,gene_start)) %>%
  GRanges()



meth.ovls <- FindOvls(regs,hg002_meth) %>%
  group_by(region,quantile,name) %>%
  summarise(meth_freq=mean(methylated_frequency)) %>%
  mutate(motif="CG")

gc.ovls <- FindOvls(regs,hg002_gc) %>%
  group_by(region,quantile,name) %>%
  summarise(meth_freq=mean(methylated_frequency)) %>%
  mutate(motif="GC")

cg.ovls <- rbind(meth.ovls,gc.ovls) %>%
  filter(motif=="CG") %>%
  spread(motif,meth_freq) %>%
  spread(region,CG) %>%
  dplyr::rename("cpg_cgi"=cgi, "cpg_gb"=gb)

gc.ovls <- rbind(meth.ovls,gc.ovls) %>%
  filter(motif=="GC") %>%
  spread(motif,meth_freq) %>%
  spread(region,GC) %>%
  dplyr::rename("gpc_cgi"=cgi, "gpc_gb"=gb)

all.ovls <- merge(cg.ovls,gc.ovls,by=c("quantile","name"))

all.ovls %>%
ggplot(aes(x=gpc_gb,y=cpg_gb,color=as.factor(quantile)))+geom_point(alpha=.5)+xlim(.2,.4)+facet_wrap(~as.factor(quantile))


my_comparisons <- list( c("1", "2"), c("1", "3"), c("1", "4"), c("1", "5"),c("2", "3"),c("2", "4"),c("2", "5"),c("3", "4"),c("3", "5"),c("4", "5"))
all.ovls %>%
 # filter(quantile %in% c(1,4)) %>%
ggplot(aes(y=gpc_gb,x=as.factor(quantile)))+geom_violin()+geom_boxplot(width=.1,outlier.shape = NA)+ stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)


ggsave(
  paste0(figs, "/HG002_violin_highvslow_GpC.pdf"),
  plot = last_plot(),
  width = 8,
  height = 5
)


all.ovls %>%
  #filter(quantile %in% c(1,4)) %>%
ggplot(aes(y=cpg_gb,x=as.factor(quantile)))+geom_violin()+geom_boxplot(width=.1,outlier.shape = NA)+ stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)

ggsave(
  paste0(figs, "/HG002_violin_highvslow_CpG.pdf"),
  plot = last_plot(),
  width = 8,
  height = 5
)

cg.bar <- all.ovls %>%
  ungroup() %>%
  group_by(quantile) %>%
  na.omit() %>%
  summarize(mean_cpg = median(as.numeric(cpg_gb)),mean_gpc = median(as.numeric(gpc_gb)),mean_prom_gpc = median(as.numeric(gpc_cgi)),mean_prom_cpg = median(as.numeric(cpg_cgi)))


chr="chrX"
library(BSgenome.HG002.chrX)
chrx.gr <- GRanges(seqinfo(BSgenome.HG002.chrX))
blocks <- genomeBlocks(BSgenome.HG002.chrX, chrs = seqnames(BSgenome.HG002.chrX), width = 15000)


rm <- read_tsv(paste0(dat, "/HG002/annotations/HG002_v0.9.orf_only.bed"), col_names = c("chr", "start", "end", "name", "len","direction")) %>%
  mutate(end=start+1) %>%
  filter(chr=="chrX") %>%
  mutate(ID = row_number()) %>%
  mutate(num=1) %>%
  GRanges()

score1 <- coverage(rm, weight="num")

binned_genes <- binnedSum(blocks, numvar = score1, "num") %>%
  as.data.frame()

genes <- ggplot(binned_genes, aes(x = start/1e6, y =num ))+theme_classic()+geom_bar(stat = "identity", position = "dodge",color="black")

ggsave(
  paste0(figs, "/ChrX_gene_density.pdf"),
  plot = genes,
  width = 9,
  height = 4
)

