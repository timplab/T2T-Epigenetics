# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/ilee_plot_utils.R")
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/methylation_R_utils.R")


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_Allchr_10kbBinned_ALL.bed") %>%
 GRanges()

hor = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621.bed"), col_names = F) %>%
  dplyr::filter(grepl("hor", X4))%>%
  mutate(status = ifelse(grepl("L", X4), "live", "dead")) %>%
  rename(X1="chr", X2="start", X3="end",X5="name") %>%
  filter(status=="live") %>%
  GRanges()

ovls <- FindOvls(all.dat,hor)

#https://stackoverflow.com/questions/6356886/detecting-dips-in-a-2d-plot
FindLowRegion <- function(x,n=length(x)/4,tol=length(x)/20,p=.9){
  nx <- length(x)
  n <-  2*(n %/% 2) + 1
  # smooth out based on means
  sx <- rowMeans(embed(c(rep(NA,n/2),x,rep(NA,n/2)),n),na.rm=T)
  # find which series are far from the mean
  rlesx <- rle((sx-x)>0)
  # construct start and end of regions
  int <- embed(cumsum(c(1,rlesx$lengths)),2)
  # which regions fulfill requirements
  id <- rlesx$value & rlesx$length > tol
  # Cut regions to be in general smaller than median
  regions <-
    apply(int[id,],1,function(i){
      i <- min(i):max(i)
      tmp <- x[i]
      id <- which(tmp < quantile(tmp,p))
      id <- min(id):max(id)
      i[id]            
    })
  # return
  unlist(regions)
}

chr1=ovls %>%
  filter(seqnames=="chr4")
Lows <- FindLowRegion(chr1$smooth)



meth <- ggplot(chr1, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)


newx <- seq_along(chr1$smooth)
newy <- ifelse(newx %in% Lows,chr1$smooth,NA)
plot(chr1$smooth, col="blue", type="l", lwd=2)
lines(newx,newy,col="red",lwd="3")
