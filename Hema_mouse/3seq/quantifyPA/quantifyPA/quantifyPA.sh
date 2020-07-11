#bulid bowtie2 index file
bowtie2-build ../../../Reference/Mm/Mus_musculus.GRCm38.dna.primary_assembly.fa GRCm38
#align
for i in $(sed '1d' ../../sample.info.txt| cut -f1 )
do
bowtie2 -x GRCm38 -U /mnt/raid63/HSC/mouse/bulk/2014science_LaraAstiasoD_chromatin/fastq/$i.fastq.gz --local  -L 25 |samtools view -bS - |samtools sort - -o ./bam/$i.sorted.bam
done

#index bamfile
 for i in $(ls *)
do
samtools index $i
done

#conbine bam files
samtools merge align/mergedbam/all.merged.bam align/bam/*sorted.bam

#separate the strand,0 pos,16 minus, 4 unmapped,the default -F is 0,unmapped reads would not appear after convert to bed 
samtools view -F 16 all.merged.bam -b |bamToBed -i /dev/stdin | cat - |sort -k 1,1 -k 2,2n  > all.pos.bed
samtools view -f 16 all.merged.bam -b |bamToBed -i /dev/stdin | cat - |sort -k 1,1 -k 2,2n  > all.neg.bed

#F seq would consider the orientation ,and we need find the 3' site of peaks
sed -i 's/+/-/' all.pos.bed
sed -i 's/-/+/' all.neg.bed
#call peaks
fseq -f 1 -l 30 -of bed ./align/mergedbam/all.pos.bed -o ./callpeaks/pos
fseq -f 1 -l 30 -of bed ./align/mergedbam/all.neg.bed -o ./callpeaks/neg

#give the strand information to call PA reads
cat ./neg/*|awk '{print$0"\t""-"}' > all.neg.peaks.bed

cat ./pos/*|awk '{print$0"\t""+"}' > all.pos.peaks.bed
cat all.neg.peaks.bed all.pos.peaks.bed |bedtools sort -i - > all.peaks.bed
grep -v "GL\|MT\|JH" all.peaks.bed > all.peaks.nomt.bed

#filter peaks
python peak_filter.py /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.dna.primary_assembly.fa ./callpeaks/all.peaks.nomt.bed /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/PAS/annoPA.sorted.bed  ./callpeaks/PA.filtered.bed

##call pa
#global PA sites
python pabulks.py ./align/mergedbam/all.merged.bam ./callpeaks/PA.filtered.bed ./all.pa/all.PA
#call pa of each sample
for i in $(sed '1d' ../sample.info.txt| cut -f1 )
do
python pabulks.py ./align/bam/$i.sorted.bam ./callpeaks/PA.filtered.bed ./all.pa/$i
done

#combine PA
for i in $(sed '1d' ../../../sample.info.txt| cut -f1 )
do
cat ../$i/* > $i.all.txt
done

cat * > ../all.PA.combined.txt
#group PA


