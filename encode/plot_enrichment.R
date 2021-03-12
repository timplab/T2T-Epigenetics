#!/usr/bin/env Rscript

### This script takes output of encode reads aligned per repeat array and generates boxplots for log2 fold enrichment
#load libraries
library(tidyverse)

# load data
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
en <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/210308_alignments/bt2.chm13v1_array_enrichmentsv2.tsv")


#set colors
censatColors =c("(CATTC)n" = "#E87C71",
                "(GAATC)n"="#E28455",
                "HOR"="#D78C32",
                "BSAT"="#E370AB",
                "CER" = "#CE9334",
                "HSAT2"="#C19935",
                "HSAT1"="#A2A638",
                "HSAT3"="#8CAC3E",
                "Low_complexity"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "MST"="#53BB73",
                "GSAT"="#55BE8D",
                "RNA"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "ACRO1"="#9400D3",
                "HSAT4"="#53B0E3",
                "SATR"="#5AA5DA",
                "CT"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "MON"="#A885BC",
                "SST"="#C378B2",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "gap-rDNA"="#ff4000",
                "TE" = "#ffbf00", 
                "TAR"= "#0080ff",
                "ACRO"="#9400D3", 
                "DHOR" = "gray")

totals <- en %>%
  filter(Array == "Total") %>%
  spread(Array, Mark) %>%
  dplyr::rename("total_control"=Control, "total_treat"=Treat, "Mark"=Total)

reads <- en %>%
  filter(Array != "Total")

sum <- merge(reads, totals, by=c("Cell type", "Mark"))

SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT")

stat.sum <- sum %>%
  group_by(`Cell type`, Mark,Array,total_control,total_treat) %>%
  summarise(array_treat=sum(Treat), array_control=sum(Control)) %>%
  summarise(enrich=log2((array_treat/array_control)*(total_control/total_treat))) %>%
  filter(Mark %in% c("H3K9me3")) %>%
  filter(Array %in% SAT)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p1 <- ggplot(stat.sum, aes(x=Mark, y=enrich))+ 
  geom_dotplot(binaxis='y', stackdir='center')+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+facet_wrap(~Array)

ggsave(
  paste0(figs, "/", "ENCODE_perCell_dotplot.pdf"),
  plot = p1,
  scale = 1,
  width = 8,
  height = 5,
)

p2 <- ggplot(stat.sum, aes(x=Array, y=enrich, fill=Array))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = 1.5)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "ENCODE_perSat_dotplot.pdf"),
  plot = p2,
  scale = 1,
  width = 8,
  height = 5,
)

p3 <- ggplot(stat.sum, aes(x=`Cell type`, y=enrich))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = .8)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")

ggsave(
  paste0(figs, "/", "ENCODE_perCellType_dotplot.pdf"),
  plot = p3,
  scale = 1,
  width = 8,
  height = 5,
)
