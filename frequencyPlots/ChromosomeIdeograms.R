require(RIdeogram)

chroms <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/Rideogram/chrom.sizes", col_names = c("Chr", "length")) %>%
  mutate(Start=0) %>%
  mutate(End=length) %>%
  select(c(Chr,Start,End))


cens <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/Rideogram/t2t-chm13.v1.0.cenSat_regions.bed", col_names = c("Chr", "Start", "End")) %>%
  mutate(Value=1)



list=c("HSAT", "BSAT", 'GSAT', "ACRO", "MON")

censat.widths <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/annotations/chm13_v1_CenSat.bed") %>%
  mutate(name = ifelse(grepl("hsat", name), "HSAT", name)) %>%
  mutate(color = ifelse(name=="HSAT", "#7570B3", name)) %>%
  mutate(name = ifelse(grepl("bsat", name), "BSAT", name)) %>%
  mutate(color = ifelse(name=="BSAT", "#66A61E", color)) %>%
  mutate(name = ifelse(grepl("GSAT", name), "GSAT", name)) %>%
  mutate(color = ifelse(name=="GSAT", "#A6761D", color)) %>%
  mutate(name = ifelse(grepl("ACRO", name), "ACRO", name)) %>%
  mutate(color = ifelse(name=="ACRO", "#666666", color)) %>%
  filter(name %in% list) %>%
  mutate(Type=name,Chr=`#chrom`,Start=chromStart,End=chromEnd, color=color, Shape="box") %>%
  select(c(Type,Shape,Chr,Start,End,color))

sst1.all <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/revision/SST1_CEN-array-mon_NONCEN-array_withLabels.txt") %>%
  mutate(Type="SST1") %>%
  group_by(`#chrom`,group,Type) %>%
  filter(group != "mon_CEN") %>%
  filter(group != "mon_nonCEN") %>%
  summarise(start=min(start), end=max(end)) %>%
  mutate(color="#1B9E77") %>%
  ungroup() %>%
  mutate(Chr=`#chrom`,Start=start,End=end, color=color, Shape="box") %>%
  select(c(Type,Shape,Chr,Start,End,color))

reps.all <- rbind(censat.widths,sst1.all)

acro <- reps.all %>%
  filter(Type=="ACRO")

ideogram(karyotype = chroms,label=acro, label_type = "marker")
convertSVG("chromosome.svg", device = "png")


acro <- reps.all %>%
  filter(Type=="HSAT")

ideogram(karyotype = chroms, label=acro, label_type = "marker")
convertSVG("chromosome.svg", device = "png")


ideogram(karyotype = chroms,overlaid = cens)
convertSVG("chromosome.svg", device = "png")
