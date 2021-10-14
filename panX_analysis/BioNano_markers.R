library(tidyverse)
library(cowplot)
sites <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.chrX.bionanoMarkers.bed") %>%
  mutate(start=start+57000000)

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


censat = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("gap", X4), "gap-rDNA", X4)) %>%
  dplyr::rename("chr" =1, "start" = 2 ,"end" =3)  %>%
    filter(chr == "chrX")


p <- ggplot(data=sites, aes(x=start/1e6, y=1))+geom_point()+theme_bw()+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(axis.text=element_text(size=18))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+labs(x="", y="")

censat.plot <- ggplot(data=censat, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
  geom_rect()+theme(legend.position="bottom") +labs(y="CenSat")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank())+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(57000000/1e6,62000000/1e6))+theme(legend.title=element_blank(),panel.border=element_blank())+theme(axis.text=element_text(size=20))

top_row <- plot_grid(p,censat.plot,ncol = 1, align="v",rel_heights = c(1/2,1/10))
