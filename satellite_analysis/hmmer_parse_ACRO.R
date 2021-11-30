dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


regs_fwd <- read_tblout(paste0(dat, "/revision_analysis/satellites/ACRO/annotate_units/ACRO_hmmsearch.tsv")) %>%
  separate(query_name, c("chr", "coords"), (":")) %>%
  separate(coords, c("num", "end"), ("-")) %>%
  rename("Orientation"=domain_accession) %>%
  select(chr, num, end,sequence_bias,best_domain_evalue,Orientation) %>%
  group_by(chr,num,end) %>%
  mutate(start=as.numeric(sequence_bias)+as.numeric(num), end=best_domain_evalue+as.numeric(num)) %>%
  mutate(len=end-start)  %>%
  ungroup() %>%
  dplyr::select(chr, start, end, len,Orientation) %>%
  mutate(Orientation=ifelse(len < 0 , "-", "+")) %>%
  mutate(len=ifelse(Orientation=="-", len*-1,len)) %>%
  mutate(ID=row_number())

regs_RC <- read_tblout(paste0(dat, "/revision_analysis/satellites/ACRO/annotate_units/ACRO_hmmsearchRC.tsv")) %>%
  separate(query_name, c("chr", "coords"), (":")) %>%
  separate(coords, c("num", "end"), ("-")) %>%
  rename("Orientation"=domain_accession) %>%
  select(chr, num, end,sequence_bias,best_domain_evalue,Orientation) %>%
  group_by(chr,num,end) %>%
  mutate(start=as.numeric(sequence_bias)+as.numeric(num), end=best_domain_evalue+as.numeric(num)) %>%
  mutate(len=end-start)  %>%
  ungroup() %>%
  select(chr, start, end, len,Orientation) %>%
  mutate(Orientation=ifelse(len < 0 , "-", "+")) %>%
  mutate(len=ifelse(Orientation=="-", len*-1,len)) %>%
  mutate(reg_start=start,reg_end=end) %>%
  mutate(start=ifelse(Orientation=="+",reg_start,reg_end)) %>%
  mutate(end=ifelse(Orientation=="+",reg_end,reg_start)) %>%
  mutate(ID=row_number()) 

write.table(regs_RC, paste0(dat, "/revision_analysis/satellites/ACRO/annotate_units/ACRO_subunits.tsv"), sep="\t", quote=F, col.names = T, row.names = F)
