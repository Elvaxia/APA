#!/bin/bash

position_list=$(cat region)

for i in $position_list
	do
		echo $i 
		sashimiplot junc --gtf /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz \
		--bam BAM.tsv \
		--junc $(echo $i|cut -f1 -d.) \
		--pa $(echo $i|cut -f3 -d.) \
		--fileout ./$(echo $i|cut -f1 -d.) \
		#--log log10\
		#--sj 10\
		--ps FR
	done

exit 0
