#working directory:/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/align
#bulid bowtie2 index file
bowtie2-build ../../../Reference/Mm/Mus_musculus.GRCm38.dna.primary_assembly.fa GRCm38

#bowtie2 alignmen
for i in $(sed '1d' ../sample.info.txt| cut -f1 )
do
bowtie2 -x GRCm38 -U /mnt/raid63/HSC/mouse/bulk/2014science_LaraAstiasoD_chromatin/fastq/$i.fastq.gz  | samtools view -bS - > $i.bam
done

#sort bam file
for i in $(sed '1d' ../sample.info.txt| cut -f1)
do
samtools sort $i.bam > $i.sorted.bam
done

#index bam file
for i in $(ls *sorted.bam)
do
samtools index $i
done
