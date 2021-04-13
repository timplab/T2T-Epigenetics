

source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)
library(BSgenome.t2t.v1.0.release)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

cen <- read_tsv(paste0(dat, "/annotations/t2t_cenRegions.v2.021621.bed"), col_names = c("chr", "start", "end","name")) %>%
  group_by(chr) %>%
  summarise(start=min(start), end=max(end)) %>%
  filter(chromosome != "chrX")


hg002_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_CpGmethylationFrequency_20kb.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrY")

hg002_flagged_bins <- read_tsv(paste0(dat, "/chm13_final_beds/HG002_flagged_CpG_CoverageBins.bed"))
chm13_flagged_bins <- read_tsv(paste0(dat, "/chm13_final_beds/CHM13_flagged_CpG_CoverageBins.bed"))

for (i in (1:length(cen$chr))){
chrom=cen$chr[i]
rstart=cen$start[i]
rend=cen$end[i]

all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_Allchr_10kbBinned_ALL.bed") %>%
  filter(seqnames == chrom) %>%
  filter(start >= rstart) %>%
  filter(end <= rend) %>%
  select(seqnames, start, end, smooth) %>%
  rename(smooth="chm13_meth")
hg002_dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/HG002_autosomes_AllCen_10kbBinned_ALL.bed") %>%
  filter(seqnames == chrom) %>%
  filter(start >= rstart) %>%
  filter(end <= rend) %>%
  select(seqnames, start, end, smooth) %>%
  rename(smooth="hg002_meth")

all <- merge(all.dat, hg002_dat, by=c("seqnames", "start","end"))
SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "gap-rDNA")

censat = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("gap", X4), "gap-rDNA", X4)) %>%
  dplyr::filter(X4 %in% SAT) %>%
  dplyr::rename("chr" =1, "start" = 2 ,"end" =3) 

meth <- ggplot()+geom_line(data=all, aes(x = (start/1e6), y= chm13_meth), color = "dodgerblue", size=1)+geom_line(data=all, aes(x = (start/1e6), y= hg002_meth), color = "purple", size=1)+ labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)

censat_sub <- censat %>%
  dplyr::filter(chr == chrom) 

censat.plot <- ggplot(data=censat_sub, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
  geom_rect()+theme(legend.position="none") +labs(y="CenSat")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = censatColors, drop = FALSE)+theme(legend.title=element_blank(),panel.border=element_blank())


chm13.flags <- chm13_flagged_bins %>%
  filter(seqnames == chrom) 


hg002.flags <- hg002_flagged_bins %>%
  filter(seqnames == chrom)

df <- data.frame(ID = c(1, 2, 3, 4, 5),
                  x = c('a', 'b', 'c', 'd', 'e'),
                  y = c(1, 1, 0, 0, 1))

if (nrow(chm13.flags)>0){
  flags.plot1 <-  ggplot()+geom_point(data=chm13.flags, aes(x = (start/1e6), y= .5), color = "dodgerblue",alpha=.5, size=1)+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(legend.title=element_blank(),panel.border=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
}  else{ 
    flags.plot1 <- ggplot(df,aes(x,y))+geom_blank()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(legend.title=element_blank(),panel.border=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }


if (nrow(hg002.flags)>0){
  flags.plot2 <-  ggplot()+geom_point(data=hg002.flags, aes(x = (start/1e6), y= .5), color = "purple",alpha=.5, size=1)+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(legend.title=element_blank(),panel.border=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
}  else{ 
  flags.plot2 <- ggplot(df,aes(x,y))+geom_blank()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(legend.title=element_blank(),panel.border=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

top_row <- plot_grid(meth,flags.plot1,flags.plot2,censat.plot,ncol = 1, align="v",rel_heights = c(1/3,1/25,1/25,1/20))
top_row
ggsave(
  paste0(figs, "/methyl_profiles/", chrom,"_PanelPlotv2maplot.pdf"),
  plot = top_row,
  scale = 1,
  width = 10,
  height = 3
)
}
