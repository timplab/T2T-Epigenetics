library(BSgenome.t2t.v1)
library(tidyverse)
library(karyoploteR)
library(liftOver)
library(ggnewscale)
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

indir="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks/cat"
end="unique.bed"
activating <- read_tsv(paste0(indir, "/activating_", end),col_names = c("chr", "start", "end", "samp", "y", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  mutate(y=1) %>%
  mutate(mark="activating")
repressing <- read_tsv(paste0(indir, "/repressing_", end),col_names = c("chr", "start", "end", "samp", "y", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  mutate(y=1) %>%
  mutate(mark="repressing")
CTCF <- read_tsv(paste0(indir, "/CTCF_", end),col_names = c("chr", "start", "end", "samp", "y", "strand", "signalValue", "pvalue", "qvalue", "peak")) %>%
  mutate(y=1) %>%
  mutate(mark="CTCF")
marks.all <- rbind(activating, repressing,CTCF)

chrom.sizes <- read_tsv(paste0("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chrom.sizes"), col_names = F) %>%
  filter(X1 != "chrM") %>%
  dplyr::rename("chromosome"=X1, "size"=X2) 


chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom.sizes[["chromosome"]] <- factor(x = chrom.sizes[["chromosome"]], 
                                      levels = chrom_order)

marks.all[["chromosome"]] <- factor(x = marks.all[["seqnames"]], 
                                 levels = chrom_order)
plot2 <- ggplot(data = chrom.sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_key)) 

plot3 <- plot2 + geom_rect(data = as.data.frame(marks.all), aes(xmin = as.numeric(chromosome) - 0.2, 
                               xmax = as.numeric(chromosome) + 0.2, 
                               ymax = start, ymin = end, fill = mark),alpha=.05)+scale_fill_manual(values=c("darkgreen", "orange", "darkred"))

ggsave(
  paste0(figs, "/overlayalpha.05_Ideogram.pdf"),
  plot = plot3,
  scale = 1,
  width = 10,
  height = 10,
)
