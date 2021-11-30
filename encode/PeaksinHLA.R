
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/figs"

all_marks <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all_HLA_intersections/HLA_numpeaks.txt", delim = "\t", col_names = F) %>%
  dplyr::select(c(X4,X5)) %>%
  dplyr::rename("num" = X4, "sample"=X5)%>%
  mutate_if(is.character, str_trim) %>%
  separate(col = sample,into = c("NBPF", "intersect","cell_line","mark"),sep = "_") %>%
  mutate_if(is.character, str_trim) %>%
  separate(col = mark,into = c("mark"),sep = "\\.") %>%
 # filter(mark %in% c('H3K4me1','H3K36me3','H3K27me3')) %>%
  mutate(num=as.numeric(num)) 

types <- read_delim("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/cell_lines.txt", delim = "\t")

all_marks <- merge(all_marks, types, by = "cell_line")

#all_marks$mark = factor(all_marks$mark, levels=c('H3K4me1','H3K36me3','H3K27me3'))

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


ggplot(data=all_marks, aes(y=num, x=cell_line, fill=cell_line))+geom_bar(stat="identity", alpha=.7,position = "identity")+ theme(text = element_text(size=20))+theme_classic()+facet_grid(~factor(mark), scales="free_x",space="free")+labs(x="Cell Line", y = "Number of MACs Peaks")+ geom_text(aes(label=num), position=position_dodge(width=0.9), vjust=-0.25)


ggsave(filename = paste0(figs, "/ENCODE_HLA_allPeaksBar.pdf"), 
       plot = last_plot(),
       width = 12,
       height = 5)


all_marks %>%
  filter(cell_line %in% c("RWPE1","RWPE2"))%>%
  filter(mark != "H3K36me3") %>%
  ggplot(aes(y=num, x=cell_line, fill=cell_line))+geom_bar(stat="identity", alpha=.7,position = "identity")+ theme(text = element_text(size=20))+theme_classic()+facet_grid(~factor(mark), scales="free_x",space="free")+labs(x="Cell Line", y = "Number of MACs Peaks")+ geom_text(aes(label=num), position=position_dodge(width=0.9), vjust=-0.25)

ggsave(filename = paste0(figs, "/ENCODE_HLA_RWPE1vs2.pdf"), 
       plot = last_plot(),
       width = 8,
       height = 5)
