#!/usr/bin/env Rscript

# install from source
#install.packages("/home/gmoney/gmoney_R_packages/makeCGI_1.3.4.tar.gz", repos = NULL, type="source")
library("makeCGI")
library(BSgenome)
library(BSgenome.t2t.v1.0.release)
library(Biostrings)
# set parameters 
.CGIoptions=CGIoptions()
.CGIoptions$rawdat.type="BSgenome"
.CGIoptions$package="BSgenome.t2t.v1.0.release"
.CGIoptions$species="t2t"
makeCGI(.CGIoptions)
