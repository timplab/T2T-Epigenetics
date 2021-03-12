#!/usr/bin Rscript

# this script takes output files generated from CG_count_chm13.R and plots bar plots for number of CG per repeat type

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

hg38 <- read_tsv("/dilithium/Data/Nanopore/Analysis/quinn/SIRV/t2t/hg38_tsv/new_hg38_counts.tsv") %>%
  mutate(asm="hg38")

chm13 <- read_tsv("/dilithium/Data/Nanopore/Analysis/quinn/SIRV/t2t/hg38_tsv/020321_chm13_counts.tsv") %>%
  mutate(asm="chm13")

repeatColors =c("DNA"="#C19936",
                "DNA?"="#C19935",
                "LINE"="#FFA500",
                "Low_complexity"="#75B043",
                "LTR"="#51B756",
                "RC"="#53BB74",
                "Retroposon"="#55BE9D",
                "RNA"="#ff4000",
                "rRNA"="#52BEBB",
                "scRNA"="#6B9AD3",
                "Simple_repeat"="#8992C9",
                "SINE"="pink",
                "snRNA"="#A886BC",
                "srpRNA"="#B67EB6",
                "Unspecified"="#C378B2",
                "tRNA"="#006400",
                "Unknown"="#BA55D3",
                "Satellite"="#53B0E4", 
                "nonrepetitive"="dodgerblue")

bind <- rbind(hg38, chm13)%>%
  mutate(types=if_else(counts_chm13 < 50000, "other", types)) %>%
  group_by(types,asm) %>%
  summarise(counts_chm13=sum(counts_chm13))


ggplot(bind, aes(x=asm, y=counts_chm13/1e6, fill= fct_reorder(bind$types, bind$counts_chm13)))+geom_bar(stat="identity", color="black")+scale_fill_manual("Repeat", values=repeatColors, drop=F)+labs(x="Assembly", y="Number CpG sites (Millions)")+theme(text = element_text(size=20))+theme_classic()

ggsave(
  paste0(figs,"/CpG_numberPerRepeat_bar.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 7
)

stat <- bind %>%
  group_by(asm) %>%
  mutate(total=sum(counts_chm13)) %>%
  group_by(asm, types) %>%
  mutate(percent=(counts_chm13/total)*100)

bp<- ggplot(stat, aes(x="", y=percent, fill=types))+
  geom_bar(width = 1, stat = "identity",color="black") +facet_wrap(~asm)
pie <- bp + coord_polar("y", start=0)+scale_fill_manual("Repeat", values=repeatColors, drop=F)+theme_void()
pie

ggsave(
  paste0(figs,"/CpG_numberPerRepeat_pie.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 7
)

stat.sum <- stat %>%
  select(-c(total, percent)) %>%
  spread(asm, counts_chm13) %>%
  mutate(dif=chm13-hg38) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(total=sum(dif)) %>%
  mutate(percent = (dif/total)*100) 

bp<- ggplot(stat.sum, aes(x="", y=percent, fill= fct_reorder(stat.sum$types, stat.sum$percent)))+
  geom_bar(width = 1, stat = "identity",color="black")
pie <- bp + coord_polar("y", start=0)+scale_fill_manual("Repeat", values=repeatColors, drop=F)+theme_void()
pie

ggsave(
  paste0(figs,"/CpG_Difference_pie.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 5,
  height = 7
)
