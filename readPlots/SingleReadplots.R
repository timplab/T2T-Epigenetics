#!/usr/bin/Rscript

# written for version 4 

# script to generate single read plots 

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(source("/home/isac/Code/ilee/plot/ilee_plot_utils.R"))
suppressPackageStartupMessages(library("ggsci"))
suppressPackageStartupMessages(source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R"))

# load data
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


# create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--region", type="character", 
                    help="Genomic region chrx:xx-xx")
parser$add_argument("-l", "--readlength", type="integer", default=50000, 
                    help="Minimum read length to include [default %(default)s]")
parser$add_argument("-m", "--maxgap", type="integer", default=20, 
                    help="Maxgap sisze for merging CpGs when plotting [default %(default)s]")
parser$add_argument("-w", "--window", type="integer", default=500, 
                    help="Number of windows for methylation frequency panel [default %(default)s]")

args <- parser$parse_args()
coords= str_split_fixed(args$region, ":",2)[1,2]
chrom=str_split_fixed(args$region, ":",2)[1,1]
rstart=str_split_fixed(coords, "-",2)[1,1]
rend=str_split_fixed(coords, "-",2)[1,2]
window=args$window
minlen=args$readlength
gap=args$maxgap

# in this dir is tabix methyl beds of all the centromeric regions for each chromosome. This command loads the per read CG calls for the centromere specified

reads <- tabix_mbed(paste0(dat, "/censat/", chrom, "_cen.mbed"),extcol = "motif",by = "read") 


size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= minlen) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

# filter for the region of interest
reg = cgcalls %>%
  dplyr::filter(start >= rstart) %>%
  dplyr::filter(end <= rend)

# for aesthetics, merge CG calls if within gap range
runs <- getRuns(reg, maxGap = gap)
# order reads based on start position
cpg_runs.ordered <- order_reads(runs)

cpg_runs_reg <- cpg_runs.ordered$x %>%
  mutate(m = ifelse(values == 1, "Methylated","Unmethylated")) %>%
  mutate(mod = "CpG")

# choose colors
pal <- pal_npg("nrc")(10)
meth_pal <- c(pal[8],pal[4]) 

# calculate methylation frequency for top panel
bis <- reg %>% 
  group_by(chrom, start) %>%
  summarise(meth_freq = mean(mcall), cov = n()) 

# bin and take rolling average to smooth methylation frequency
bis$cut = cut(bis$start, breaks=window)

meth.stat <- bis %>%
  group_by(cut) %>%
  summarise(med = median(meth_freq), top = quantile(meth_freq, 0.75), bot = quantile(meth_freq, 0.25), n_genes = length(meth_freq)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  mutate(med_smooth=rollmean(med, 10, NA))

#plot
meth <- ggplot(meth.stat, aes(x = min, y= med_smooth))+geom_line()+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), panel.spacing = unit(2, "lines"))


g <- ggplot(cpg_runs_reg,aes(xmin = start/1e6, xmax = end/1e6, ymin = ymin, ymax = ymax)) + geom_rect(data = cpg_runs.ordered$bounds, fill = "grey80") + 
  geom_rect(data=cpg_runs_reg,aes(fill = m))  +
  #    geom_vline(xintercept = 127638255, linetype == "dashed") +
  scale_fill_manual(name = "State", values = meth_pal) + theme(axis.text.y = element_blank(),  legend.position= "bottom",axis.ticks.y = element_blank(), panel.spacing = unit(2, "lines"))

top_row <- plot_grid(meth,g, align="v", ncol=1,rel_heights = c(1/10, 1/2))


ggsave(paste0(chrom,  "-", rstart, "-", rend, "_singleRead.pdf"),
  plot = top_row,
  width = 8,
  height = 8)

meth.plot <- meth.stat %>%
  select(c(min,med_smooth)) %>%
  column_to_rownames("min") %>%
  na.omit() %>%
  as.matrix()

colnames(meth.plot) <- NULL

pdf(paste0(chrom,  "-", rstart, "-", rend,"_autocorrelationPlot.pdf")) 
acf(meth.plot,lag=100)
dev.off()
