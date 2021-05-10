#!/usr/bin/Rscript

# Creator: Paul W. Hook
# Date: March 8, 2021

# load libraries
library(SRAdb)
library(DBI)

# grab SRA database
srafile <- getSRAdbFile()

# Connect to database
con=dbConnect(RSQLite::SQLite(),srafile)


listSRAfile('SRP028804',con)
sralist <- listSRAfile('SRP028804',con)
sralist

# Write out data
write.table(x=sralist, 
file="/mithril/Data/NGS/projects/rrbs_t2t/GSE49828_sra_list.txt",
quote=FALSE,col.names=TRUE,
row.names=FALSE,sep="\t")

# Grab the SRR files names
srr <- sralist$run
srr
srr <- as.data.frame(sralist$run)
write.table(x=srr,
file="/mithril/Data/NGS/projects/rrbs_t2t/GSE49828_srr_list.txt",
quote=FALSE,
col.names=FALSE,
row.names=FALSE,
sep="\t")
