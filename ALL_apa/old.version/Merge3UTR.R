#merge 3UTR of gene from gencode gtf
library(data.table)
library(dplyr)
all.UTR <- fread('/mnt/data1/xiaoxia/ALL_apa/UTR3.ref.txt', header = T)

#merge all transcripts 3UTR to generate one 3'UTR
setkey(all.UTR, gene_id)

left <- all.UTR[,lapply(.SD, min), by = gene_id, .SDcol = 'start']
right <- all.UTR[,lapply(.SD, max), by = gene_id, .SDcol = 'end']
strand <- all.UTR[,lapply(.SD, unique), by = gene_id, .SDcol = 'strand']
chr <- all.UTR[,lapply(.SD, unique), by = gene_id, .SDcol = 'chr']
gene.name <- all.UTR[,lapply(.SD, unique), by = gene_id, .SDcol = 'gene_name']


utr3.bed <- data.frame(
  chr = chr$chr,
  start = left$start,
  end = right$end,
  id = paste(chr$gene_id, gene.name$gene_name,chr$chr,strand$strand,sep = '|'),
  sc = 0,
  strand = strand$strand
)

##TODO :::extended 3UTR 10kb


write.table(utr3.bed, file = 'GRCh38.3utr.unextended.bed', sep = '\t',
            col.names = F, row.names = F, quote = F)
