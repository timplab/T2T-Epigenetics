library(tidyverse)
library(ggrepel)

read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/H3K9me3_intersections/all_peaks.txt", col_names = c("peaks", "cell_line", "region")) %>%
  group_by(cell_line) %>%
  mutate(total=sum(peaks)) %>%
  group_by(cell_line, region) %>%
  mutate(frac=peaks/total) %>%
  ggplot(aes(x="", y=frac, fill=region)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+facet_wrap(~cell_line)+theme_void()+
  geom_label_repel(aes(label = peaks), size=5, show.legend = F, nudge_x = 1) +
  guides(fill = guide_legend(title = "Group"))

ggsave(
  "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/H3K9me3_intersections/H3K9me3_pie.pdf",
  plot = last_plot(),
  width = 9,
  height = 5
)

