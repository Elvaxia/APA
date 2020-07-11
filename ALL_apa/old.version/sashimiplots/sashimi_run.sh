#!/bin/bash

position_list=$(cat test.region)

for i in $position_list
	do
		echo $i 
		sashimiplot junc --gtf /mnt/raid61/APA/Project/ALL_apa/ref/gtf.sorted.gz \
		--bam BAM.tsv \
		--junc $(echo $i|cut -f1 -d.) \
		--pa $(echo $i|cut -f3 -d.) \
		--fileout ./TEST/$(echo $i|cut -f1 -d.) \
		#--log log10\
		#--sj 10\
		#--ps FR
	done

exit 0
