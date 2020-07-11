#filter MAPQ >=20
for i in $(sed '1d' sample.info.txt| cut -f1 )
do
samtools view -Bh -q 20 ./align/$i.sorted.bam > ./align_filtered/$i.filtered.bam
done
#sort filtered bam files
for i in $(sed '1d' sample.info.txt| cut -f1 )
do
samtools sort ./align_filtered/$i.filtered.bam > ./align_filtered/$i.filtered.sorted.bam
done
#index filtered bam files
for i in $(sed '1d' sample.info.txt| cut -f1 )
do
samtools index ./align_filtered/$i.filtered.sorted.bam
done

# rm *filtered.bam
#generate bed file

for i in $(sed '1d' ../sample.info.txt| cut -f1 )
do
bamToBed -i $i.filtered.sorted.bam |bedtools sort -i - > ./bed/$i.filtered.sorted.bed
done
