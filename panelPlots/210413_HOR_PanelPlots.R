# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/ilee_plot_utils.R")
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/methylation_R_utils.R")


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
list=c("chr3","chr5", "chr19")

# chr19 coords 24500000 31000000
# chr5 coords 
# chr 3 coords 90000000 97000000
cen = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("gap", X4), "gap-rDNA", X4)) %>%
  dplyr::filter(X4 == "HOR") %>%
  dplyr::rename("chrom" =1, "start" = 2 ,"end" =3) %>%
  group_by(chrom) %>%
  summarize(start=min(start), end=max(end))

for (i in (1:length(cen$chrom))){
  chrom=cen$chrom[i]
  rstart=cen$start[i]-100000
  rend=cen$end[i]+100000
  
  all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_Allchr_10kbBinned_ALL.bed") %>%
    filter(seqnames == chrom) %>%
    filter(start >= rstart) %>%
    filter(end <= rend)
  
  
  repeatColors =c("DNA"="#C19936",
                  "DNA?"="#C19935",
                  "LINE"="#FFA500",
                  "Low_complexity"="#75B043",
                  "LTR"="#51B756",
                  "RC"="#53BB74",
                  "Retroposon"="#55BE9D",
                  "RNA"="#ff4000",
                  "rRNA"="#52BEBB",
                  "scRNA"="#6B9AD3",
                  "Simple_repeat"="#8992C9",
                  "SINE"="#9A8AC1",
                  "snRNA"="#A886BC",
                  "srpRNA"="#B67EB6",
                  "Unspecified"="#C378B2",
                  "tRNA"="#006400",
                  "Unknown"="#BA55D3",
                  "Satellite"="#53B0E4")
  
  defaultColor = "#000080"
  
  censatColors =c("TE" = "#E87C71",
                  "MER"="#E28455",
                  "HOR"="#D78C32",
                  "BSAT"="#E370AB",
                  "CER" = "#CE9334",
                  "HSAT2"="#C19935",
                  "HSAT1"="#A2A638",
                  "HSAT3"="#8CAC3E",
                  "L1"="#75B042",
                  "LSAU"="#54B346",
                  "LTR"="#51B756",
                  "MST"="#53BB73",
                  "GSAT"="#55BE8D",
                  "GSATII"="#54C0A5",
                  "rRNA"="#52BEBB",
                  "SAR"="#51BDCE",
                  "novel"="#9400D3",
                  "HSAT4"="#53B0E4",
                  "SATR"="#5AA5DA",
                  "CT"="#6B9AD2",
                  "HERV"="#8992C8",
                  "MSAT"="#9A8AC2",
                  "MON"="#A885BC",
                  "SST"="#C378B2",
                  "HSAT5"="#ED72A5",
                  "HSAT6"="#EF768C", 
                  "gap-rDNA"="#ff4000",
                  "L1" = "#ffbf00", 
                  "TAR"= "#0080ff",
                  "ACRO"="#9400D4",
                  "Alu"="#9A8AC3", 
                  "DHOR"="grey")
  
  SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "gap-rDNA")
  censat = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = F) %>%
    mutate(name = ifelse(grepl("gap", X4), "gap-rDNA", X4)) %>%
    dplyr::filter(X4 %in% SAT) %>%
    dplyr::rename("chr" =1, "start" = 2 ,"end" =3) 
  
  meth <- ggplot(all.dat, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)
 
   cov <- ggplot(all.dat, aes(x = (start/1e6), y= ont_cov))+geom_histogram(stat = "identity") + labs(y="Coverage")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
   
  censat_sub <- censat %>%
    dplyr::filter(chr == chrom) 
  
  censat.plot <- ggplot(data=censat_sub, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
    geom_rect()+theme(legend.position="none") +labs(y="CenSat")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = censatColors, drop = FALSE)+theme(legend.title=element_blank(),panel.border=element_blank())
  
  hor = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621.bed"), col_names = F) %>%
    dplyr::filter(grepl("hor", X4))%>%
    mutate(status = ifelse(grepl("L", X4), "live", "dead"))%>%
    filter(X1 == chrom) 
  
  hor.plot <- ggplot()+geom_rect(data=hor, aes(xmin=(X2/1e6),xmax=(X3/1e6),ymin=0,ymax=.1, fill=status, color=status, alpha=.5))+geom_rect()+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "none")+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))
  
  top_row <- plot_grid(meth,cov,hor.plot,censat.plot,ncol = 1, align="v",rel_heights = c(1/3,1/10,1/25,1/20))
  top_row
  ggsave(
    paste0(figs, "/methyl_profiles/", chrom,"_HOR_PanelPlotv2Cov.pdf"),
    plot = top_row,
    scale = 1,
    width = 10,
    height = 3
  )

}
