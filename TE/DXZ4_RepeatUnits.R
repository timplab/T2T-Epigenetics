units <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/DXZ4/DXZ4_nhmmer_parsed.tsv", 
                  col_names = c("chr", "start", "end", "strand"))

arrow <- ggplot()+ geom_segment(data = units, aes(x = start, xend = end, y = 0, yend = 0), 
                                colour = "black", 
                                arrow = arrow())+theme_void()


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/revision"
ggsave(
  paste0(figs,"/DXZ4_repeat_units.pdf"),
  plot = arrow,
  scale = 1,
  width = 10,
  height = 4,
)

