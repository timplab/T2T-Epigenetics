#!/bin/sh

# I have to manually put stuff together so here it is. I hate it
# PWH 02/26/2021

# 8-cell1
~/software/Bismark-0.22.2/bismark2bedGraph -o 8-cell1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950984_bismark_bt2_pe.txt.gz CpG_context_SRR950985_bismark_bt2_pe.txt.gz &> processed_meth_calls/eightcell1_cpg.log
# 8-cell2
~/software/Bismark-0.22.2/bismark2bedGraph -o 8-cell2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950986_bismark_bt2_pe.txt.gz CpG_context_SRR950987_bismark_bt2_pe.txt.gz &> processed_meth_calls/eightcell2_cpg.log
# 8-cell3
~/software/Bismark-0.22.2/bismark2bedGraph -o 8-cell3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950988_bismark_bt2_pe.txt.gz CpG_context_SRR950989_bismark_bt2_pe.txt.gz &> processed_meth_calls/eightcell3_cpg.log
# ICM1
~/software/Bismark-0.22.2/bismark2bedGraph -o ICM1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950990_bismark_bt2_pe.txt.gz CpG_context_SRR950991_bismark_bt2_pe.txt.gz &> processed_meth_calls/icm1_cpg.log
# Morula1
~/software/Bismark-0.22.2/bismark2bedGraph -o Morula1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950996_bismark_bt2_pe.txt.gz CpG_context_SRR950997_bismark_bt2_pe.txt.gz &> processed_meth_calls/morula1_cpg.log
# Morula2
~/software/Bismark-0.22.2/bismark2bedGraph -o Morula2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950998_bismark_bt2_pe.txt.gz CpG_context_SRR950999_bismark_bt2_pe.txt.gz &> processed_meth_calls/morula2_cpg.log
# Morula3
~/software/Bismark-0.22.2/bismark2bedGraph -o Morula3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951000_bismark_bt2_pe.txt.gz CpG_context_SRR951001_bismark_bt2_pe.txt.gz &> processed_meth_calls/morula3_cpg.log
# Postimplantation1
~/software/Bismark-0.22.2/bismark2bedGraph -o Postimplantation1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951002_bismark_bt2_pe.txt.gz CpG_context_SRR951003_bismark_bt2_pe.txt.gz &> processed_meth_calls/postimplant1_cpg.log
# Postimplantation2
~/software/Bismark-0.22.2/bismark2bedGraph -o Postimplantation2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951004_bismark_bt2_pe.txt.gz CpG_context_SRR951005_bismark_bt2_pe.txt.gz &> processed_meth_calls/postimplant2_cpg.log
# Postimplantation3
~/software/Bismark-0.22.2/bismark2bedGraph -o Postimplantation3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951006_bismark_bt2_pe.txt.gz CpG_context_SRR951007_bismark_bt2_pe.txt.gz &> processed_meth_calls/postimplant3_cpg.log
# TE1
~/software/Bismark-0.22.2/bismark2bedGraph -o TE1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951012_bismark_bt2_pe.txt.gz CpG_context_SRR951013_bismark_bt2_pe.txt.gz &> processed_meth_calls/te1_cpg.log
# TE2
~/software/Bismark-0.22.2/bismark2bedGraph -o TE2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951014_bismark_bt2_pe.txt.gz CpG_context_SRR951015_bismark_bt2_pe.txt.gz &> processed_meth_calls/te2_cpg.log

# 1st-PB1,SRR950976
~/software/Bismark-0.22.2/bismark2bedGraph -o 1st-PB1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950976_bismark_bt2_pe.txt.gz &> processed_meth_calls/1st-PB1_cpg.log

# 1st-PB2,SRR950977
~/software/Bismark-0.22.2/bismark2bedGraph -o 1st-PB2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950977_bismark_bt2_pe.txt.gz &> processed_meth_calls/1st-PB2_cpg.log
# 1st-PB3,SRR951019
~/software/Bismark-0.22.2/bismark2bedGraph -o 1st-PB3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951019_bismark_bt2_pe.txt.gz &> processed_meth_calls/1st-PB3_cpg.log
# 2-cell1,SRR950980
~/software/Bismark-0.22.2/bismark2bedGraph -o 2-cell1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950980_bismark_bt2_pe.txt.gz &> processed_meth_calls/2-cell1_cpg.log
# 2-cell2,SRR950981
~/software/Bismark-0.22.2/bismark2bedGraph -o 2-cell2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950981_bismark_bt2_pe.txt.gz &> processed_meth_calls/2-cell2_cpg.log
# 2nd-PB1,SRR950978
~/software/Bismark-0.22.2/bismark2bedGraph -o 2nd-PB1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950978_bismark_bt2_pe.txt.gz &> processed_meth_calls/2nd-PB1_cpg.log
# 2nd-PB2,SRR950979
~/software/Bismark-0.22.2/bismark2bedGraph -o 2nd-PB2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950979_bismark_bt2_pe.txt.gz &> processed_meth_calls/2nd-PB2_cpg.log
# 4-cell1,SRR950982
~/software/Bismark-0.22.2/bismark2bedGraph -o 4-cell1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950982_bismark_bt2_pe.txt.gz &> processed_meth_calls/4-cell1_cpg.log
# 4-cell2,SRR950983
~/software/Bismark-0.22.2/bismark2bedGraph -o 4-cell2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950983_bismark_bt2_pe.txt.gz &> processed_meth_calls/4-cell2_cpg.log
# ICM2,SRR950992
~/software/Bismark-0.22.2/bismark2bedGraph -o ICM2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950992_bismark_bt2_pe.txt.gz &> processed_meth_calls/ICM2_cpg.log
# ICM3,SRR950993
~/software/Bismark-0.22.2/bismark2bedGraph -o ICM3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950993_bismark_bt2_pe.txt.gz &> processed_meth_calls/ICM3_cpg.log
# MII-Oocyte1,SRR950994
~/software/Bismark-0.22.2/bismark2bedGraph -o MII-Oocyte1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950994_bismark_bt2_pe.txt.gz &> processed_meth_calls/MII-Oocyte1_cpg.log
# MII-Oocyte2,SRR950995
~/software/Bismark-0.22.2/bismark2bedGraph -o MII-Oocyte2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR950995_bismark_bt2_pe.txt.gz &> processed_meth_calls/MII-Oocyte2_cpg.log
# Sperm1.SRR951008
~/software/Bismark-0.22.2/bismark2bedGraph -o Sperm1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951008_bismark_bt2_pe.txt.gz &> processed_meth_calls/Sperm1_cpg.log
# Sperm2,SRR951009
~/software/Bismark-0.22.2/bismark2bedGraph -o Sperm2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951009_bismark_bt2_pe.txt.gz &> processed_meth_calls/Sperm2_cpg.log
# Sperm3,SRR951010
~/software/Bismark-0.22.2/bismark2bedGraph -o Sperm3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951010_bismark_bt2_pe.txt.gz &> processed_meth_calls/Sperm3_cpg.log
# Sperm4,SRR951011
~/software/Bismark-0.22.2/bismark2bedGraph -o Sperm4_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951011_bismark_bt2_pe.txt.gz &> processed_meth_calls/Sperm4_cpg.log
# TE3,SRR951016
~/software/Bismark-0.22.2/bismark2bedGraph -o TE3_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951016_bismark_bt2_pe.txt.gz &> processed_meth_calls/TE3_cpg.log
# Zygote1,SRR951017
~/software/Bismark-0.22.2/bismark2bedGraph -o Zygote1_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951017_bismark_bt2_pe.txt.gz &> processed_meth_calls/Zygote1_cpg.log
# Zygote2,SRR951018
~/software/Bismark-0.22.2/bismark2bedGraph -o Zygote2_cpg.bedgraph.gz --buffer_size 20G --remove_spaces --dir processed_meth_calls CpG_context_SRR951018_bismark_bt2_pe.txt.gz &> processed_meth_calls/Zygote2_cpg.log
