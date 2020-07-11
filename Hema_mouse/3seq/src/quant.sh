for i in $(sed '1d' sample.info.txt| cut -f1 )
do
bedtools intersect -a ./callpeaks/filter_peaks/final.filtered.all.bed -b ./align_filtered/bed/$i.filtered.sorted.bed -c -s > ./quant/$i.counts.bed
done
