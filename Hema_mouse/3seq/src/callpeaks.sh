#call peaks from conbined bam files
samtools merge merged_bam/all.merged.bam ../align_filtered/*sorted.bam
#tag 0 means align to the positive strand/16 minus strand,as samtools view default is 0, so use -F just exclude the minus 
samtools view -F 16 all.merged.bam -b |bamToBed -i /dev/stdin | cat - |sort -k 1,1 -k 2,2n  > all.pos.bed
samtools view -f 16 all.merged.bam -b |bamToBed -i /dev/stdin | cat - |sort -k 1,1 -k 2,2n  > all.neg.bed

fseq -f 1 -l 30 -of bed ./merged_bam/all.pos.bed -o global_peaks/pos
fseq -f 1 -l 30 -of bed ./merged_bam/all.neg.bed -o global_peaks/neg
cat ./neg/*|awk '{print$0"\t""-"}' > all.neg.peaks.bed

cat ./pos/*|awk '{print$0"\t""+"}' > all.pos.peaks.bed

#filter enrichment score >=1
awk '{if($5>=1) print$0}' ../global_peaks/all.pos.peaks.bed |bedtools sort -i - > filtered.score.peaks.pos.sorted.bed

awk '{if($5>=1) print$0}' ../global_peaks/all.neg.peaks.bed |bedtools sort -i - > filtered.score.peaks.neg.sorted.bed
#filter reads >= 50
bedtools intersect -a filtered.score.peaks.neg.sorted.bed -b ../merged_bam/all.neg.bed -c  |awk '{if($7>=50) print$0}' - > final.filtered.neg.bed
bedtools intersect -a filtered.score.peaks.pos.sorted.bed -b ../merged_bam/all.pos.bed -c  |awk '{if($7>=50) print$0}' - > final.filtered.pos.bed




