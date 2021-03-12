
library(tidyverse)
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
stat.sum <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/cut_and_run/H3K9me3_51merintersections_log2FC")


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

ggplot(stat.sum, aes(x=name, y=enrich, fill=name))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',dotsize = 1.5)+theme_classic()+geom_hline(yintercept = 0, linetype="dashed")+scale_fill_manual(values=censatColors)

ggsave(
  paste0(figs, "/", "CHm13H3K9Me3_enrichment.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 5,
)