#!/bin/bash 
repo=/home/gmoney/repos/rerio
guppy_path=/home/gmoney/repos/ont-guppy/bin/guppy_basecall_server
ref=/atium/Data/old_mithril_ref/human38/GRCH38.fa
#guppy_basecaller -i ../fast5_pass -s basecalled_fast5s \
#	-d ${repo}/basecall_models/ \
#	-x "cuda:0" \
#	-c res_dna_r941_min_modbases-all-context_v001.cfg

megalodon ../fast5_pass --mod-motif Z CG 0 --devices 0 --processes 40 \
	--outputs basecalls mappings per_read_mods mod_mappings\
	--write-mods-text \
	--reference $ref \
	--overwrite \
	--guppy-server-path $guppy_path \
	--verbose-read-progress 3 \
	--mod-map-base-conv C T --mod-map-base-conv Z C \
	--mod-map-emulate-bisulfite
#-guppy-params "-d /home/gmoney/repos/rerio/basecall_models/" \
#--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
