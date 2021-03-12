#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
data_fname = args[1]
fig_fname = args[2]

data <- read.table(data_fname, header=T)
data = data[data$group != "total",]
data$delta = data$chm13v1 / data$GRCh38p13

groups = c("mapped", "filtered", "dedup", "kmer")
genomes = c("GRCh38p13", "chm13v1")

colors = c("mapped"="#ED72A5",
           "filtered"="#55BE8D",
           "dedup"="#53B0E4",
           "kmer"="#D78C32")
shapes = c("mapped"="triangle",
           "filtered"="square",
           "dedup"="circle")
label_dict = c("mapped"="All",
               "filtered"="Paired primary",
               "dedup"="Unique primary",
               "kmer"="K-mer filtered")

data$color=colors[data$group]
data$label=""

labels = c()
for (i in (1:length(groups))){
  meanval = mean(data$delta[data$group == groups[i]])
  labels = c(labels, paste0(label_dict[groups[i]], " (", round(meanval, 3), ")", sep=""))
}

plot1 = ggplot(data, aes(x=GRCh38p13, y=delta, color=group)) + geom_point(size=3) + geom_line(aes(group=expID))
plot1 = plot1 + geom_hline(yintercept=1)
plot1 = plot1 + labs(x="# Mapped Reads GRCh38p13", y="% change in chm13v1")
plot1 = plot1 + theme_classic(base_size = 13)
plot1 = plot1 + scale_color_manual(labels=labels, breaks=groups, values=colors)
plot1 = plot1 + theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position="bottom")


ggsave(
  fig_fname,
  plot = plot1,
  scale = 1,
  width = 7.2,
  height = 7.2
)

