library(tidyverse)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"
all_marks <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds/all_peaks_num.txt", delim = "\t", col_names = F) %>%
  dplyr::select(c(X5,X6)) %>%
  dplyr::rename("num" = X5, "sample"=X6)%>%
  mutate_if(is.character, str_trim) %>%
  separate(col = sample,into = c("cell_line", "rest"),sep = "_") %>%
  mutate_if(is.character, str_trim) %>%
  separate(col = rest,into = c("mark", "asm"),sep = "\\.") %>%
  mutate(num=as.numeric(num)) %>%
  filter(asm != "GRCh38p13LO")

types <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/cell_lines.txt", delim = "\t")

all_marks <- merge(all_marks, types, by = "cell_line")

all_marks$mark = factor(all_marks$mark, levels=c('H3K27ac','H3K4me1','H3K36me3','H3K4me1','H3K4me3','H3K27me3','H3K9me3','CTCF'))

all_marks$cell_line = factor(all_marks$cell_line, levels=c("HL-60",
                                                           "MG63",
                                                           "SJCRH30",
                                                           "SJSA1",
                                                           "BE2C",
                                                           "brain",
                                                           "C4-2B",
                                                           "Caco-2",
                                                          "HAP-1",
                                                           "22Rv1",
                                                           "prostate",
                                                           "RWPE1",
                                                           "RWPE2",
                                                           "VCaP"))

ggplot(data=all_marks, aes(y=num, x=cell_line, fill=asm))+geom_bar(stat="identity", alpha=.7,position = "identity")+ theme(text = element_text(size=20))+theme_classic()+facet_grid(~factor(mark, levels=c('H3K27ac','H3K4me1','H3K36me3','H3K4me3','H3K27me3','H3K9me3',"CTCF")), scales="free_x",space="free")+labs(x="Cell Line", y = "Number of MACs Peaks")

ggsave(filename = paste0(figs, "/GRCh38_chm13peaks_vertical.pdf"), plot = last_plot(),
       width = 10,
       height = 3)


all_marks <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks/beds/all_peaks_num.txt", delim = "\t", col_names = F) %>%
  dplyr::select(c(X5,X6)) %>%
  dplyr::rename("num" = X5, "sample"=X6)%>%
  mutate_if(is.character, str_trim) %>%
  separate(col = sample,into = c("cell_line", "rest"),sep = "_") %>%
  mutate_if(is.character, str_trim) %>%
  separate(col = rest,into = c("mark", "asm"),sep = "\\.") %>%
  mutate(num=as.numeric(num)) 

stat.sum <- all_marks %>%
  spread(asm, num) %>%
  group_by(mark) %>%
  dplyr::summarize(chm13_all = sum(chm13v1),chm13_new = sum(chm13v1-GRCh38p13),LO_loss=sum(GRCh38p13-GRCh38p13LO) ,shared=sum(GRCh38p13LO), n=n(),perc=(chm13_new/chm13_all)*100) %>%
  distinct()


stat.sum <- all_marks %>%
  spread(asm, num) %>%
  group_by(mark) %>%
  dplyr::summarize(chm13_all = sum(chm13v1),hg38_all= sum(GRCh38p13),chm13_new = sum(chm13v1-GRCh38p13),LO_loss=sum(GRCh38p13-GRCh38p13LO), n=n(),perc=(chm13_new/chm13_all)*100) %>%
  distinct()

novel <-  read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/unique_peaks/NewPeaks.txt", delim = "\t", col_names = F) %>%
  dplyr::select(c(X5,X6)) %>%
  dplyr::rename("num" = X5, "sample"=X6)%>%
  mutate_if(is.character, str_trim) %>%
  separate(col = sample,into = c("cell_line", "rest"),sep = "_") %>%
  mutate_if(is.character, str_trim) %>%
  separate(col = rest,into = c("mark"),sep = "\\.") %>%
  mutate(num=as.numeric(num)) %>%
  group_by(mark) %>%
  summarise(total=sum(num), n=n())
