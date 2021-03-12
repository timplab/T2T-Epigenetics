setwd("/home/gmoney/gmoney_R_packages")
forgeBSgenomeDataPkg("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/seed_file.txt")

# on command line 
# /usr/bin/R CMD build BSgenome.t2t.v1.0.release
# /usr/bin/R CMD check BSgenome.t2t.v1.0.release.tar.gz
# /usr/bin/R CMD INSTALL BSgenome.t2t.v1.0.release_1.0.0.tar.gz

# check 

library(BSgenome.t2t.v1.0.release)

#woot
