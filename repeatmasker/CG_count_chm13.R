# this script takes repeatmasker output bed files and counts number of CG sites within each repeat type
# author Quinn Hauck 2020

##load t2t release as BSgenome object
load_all('~/R/BSgenome.t2t.v1.0.release/')
library(BSgenome.t2t.v1.0.release)
genome <- BSgenome.t2t.v1.0.release

##install("BSgenome.Hsapiens.UCSC.hg38")
##genome <- getBSgenome('hg38')
##change above back!!!!!!
##load repeatmasker BED file as a GRanges object
repeats <- readBed('/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations/chm13.draft_v1.0_plus38Y_repeatmasker.out.bed')

## establish pattern of interest, find total count in T2T genome
cg <- 'cg'
cg_count_chm13 <- 0
for (i in 1:22) {
  temp_chrm_name <- paste0('chr', i)
  temp_chrm <- genome[[temp_chrm_name]]
  cg_count_chm13 <<- cg_count_chm13 + countPattern(cg, temp_chrm)
}
cg_count_chm13 <- cg_count_chm13 + countPattern(cg, genome$chrX)

## Filter Granges object into dataframe of only desired variables

repeats_df <- as.data.frame(repeats)
repeats_df <- repeats_df %>%
  dplyr::select(seqnames, start, end, strand, name, thickStart)
colnames(repeats_df) <- c('chrm', 'start', 'end', 'strand', 'name', 'type')
##repeats_df <- repeats_df %>%
  ##filter(chrm == 'chrMT')

## Create a vector of all types of repeats
types <- unique(repeats_df$type)
counts_chm13 <- vector(mode = 'numeric', length = length(types))
##iterate through the dataframe, adding a CG count for every element in the repeatmasker
##first separate by chromosome
for (j in 1:23){
  if (j!=23){
    temp_chr <- paste0('chr', j)
  unique_chr_df <- repeats_df %>%
    filter(chrm == temp_chr)
  
  indi_chrm_values <- vector(mode = 'numeric', length = length(types))
  }
  else{
    temp_chr <- 'chrX'
    unique_chr_df <- repeats_df %>%
      filter(chrm == temp_chr)

    indi_chrm_values <- vector(mode = 'numeric', length = length(types))
  }
##iterate through each type of sequence
for (i in 1:length(types)){
  specific_repeat_df <- unique_chr_df %>%
    filter(type == types[i])
  
  ## add a CG count to corresponding vector entry in counts for each element of selected type
  if(nrow(specific_repeat_df != 0)){
    for (k in 1:(nrow(specific_repeat_df))){
      temp_start <- as.numeric((specific_repeat_df$start)[k])
      temp_end <- as.numeric((specific_repeat_df$end)[k])
      temp_loc <- (genome[[temp_chr]])[temp_start:temp_end]
      temp_value <- countPattern(cg, temp_loc)
      counts_chm13[i] <- counts_chm13[i] + temp_value
      indi_chrm_values[i] <- indi_chrm_values[i] + temp_value
    }
  }
}
  temp_chr_count <- countPattern(cg, genome[[temp_chr]])
  temp_types_chm13 <- append(types, 'nonrepetitive')
  indi_chrm_values <- append(indi_chrm_values, temp_chr_count-sum(indi_chrm_values))
  indi_chrm_df <- data.frame(temp_types_chm13, indi_chrm_values)
  file_name <- paste0('/dilithium/Data/Nanopore/Analysis/quinn/SIRV/t2t/chm13_tsv/', temp_chr)
  write_tsv(indi_chrm_df, file_name)

}


types <- append(types, 'nonrepetitive')
counts_chm13 <- append(counts_chm13, cg_count_chm13-sum(counts_chm13))
full_breakdown <- data.frame(types, counts_chm13)

##remove any elements that make up less than 1%, group into other category
other <- 0
threshold <- 0.01*cg_count_chm13
storage <- c()
for (i in 1:nrow(full_breakdown)){
  if(counts_chm13[i] < threshold){
    other <- other + counts_chm13[i]
    storage <- as.vector(append(storage, i))
  }
}
full_breakdown <- full_breakdown[-storage,]

other_vector <- c('other', as.numeric(other))
full_breakdown <- rbind(full_breakdown, other_vector)

full_breakdown$counts_chm13 <- as.numeric(as.character(full_breakdown$counts_chm13))
full_breakdown <- full_breakdown %>%
  mutate(value = counts_chm13/cg_count_chm13)

bar <- ggplot(full_breakdown, aes(x="", y=value, fill=types))+
  geom_bar(width = 1, stat = "identity")
pie <- bar + coord_polar("y", start=0)
pie

df <- cbind(types, counts_chm13)
write_tsv(as.data.frame(df), '/dilithium/Data/Nanopore/Analysis/quinn/SIRV/t2t/hg38_tsv/020321_chm13_counts.tsv')
