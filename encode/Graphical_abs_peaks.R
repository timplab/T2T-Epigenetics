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

sum.marks <-all_marks %>%
  group_by(mark,asm) %>%
  summarise(num=mean(num))%>%
  mutate(type=ifelse(mark %in% c('H3K4me1','H3K27ac','H3K36me3',"CTCF",'H3K4me3'), "activating", "repressing"))

p1 <- ggplot(data=sum.marks, aes(y=num, x=asm, color=mark, group=mark,shape=type))+geom_point()+geom_path(position = "identity")+ theme(text = element_text(size=20))+theme_classic()+labs( y = "Number of MACs Peaks") + scale_y_break(c(6e4, 9.9e4))+theme(legend.position = "bottom")

ggsave(filename = paste0(figs, "/GRCh38_chm13peaks_vertical_lines.pdf"), plot = p1,
       width = 5,
       height = 5)
