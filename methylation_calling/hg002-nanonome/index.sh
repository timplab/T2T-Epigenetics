#!/bin/bash 

path=/data/201201_HG002_nanonome_SREXL/
fastq=${path}/201201_HG002_nanonome_SREXL_4-5.fastq
~/repos/f5c-v0.5/f5c_x86_64_linux index -d ${path}/201201_HG002_nanonome_SREXL_4/20201201_2120_2-E3-H3_PAF27413_e7b711ca/fast5_pass -d ${path}/201201_HG002_nanonome_SREXL_5/20201201_2120_2-A9-D9_PAF25689_9f652d26/fast5_pass ${fastq} --iop 3 -t 48 
