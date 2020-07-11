#!/usr/bin/env Rscript


purpose <- 'Extract and annotate the intron/3UTR from Ensembl gtf file.'

##------------------------------------------------------
## setup Rscript
suppressPackageStartupMessages({library(optparse)})

Usage <- paste('%prog [OPTIONS] \n', purpose)
option.list <-
  list(make_option(c('-o', '--outdir'), default = '.', type = 'character',
                   help = 'The directory to output the results to'),
       make_option(c('-f', '--gtf'), default = NULL, type = 'character',
                   help = 'gtf file download from Ensembl')
  )

opt_parser = OptionParser(usage=Usage, option_list=option.list)
opt = parse_args(opt_parser)
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("gtf file must be supplied .n", call.=FALSE)
}

##------------------------------------------------------
## script start from here
suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(annotatr)
  library(dplyr)
})


gtf_file <- opt$gtf


#gtf_file <- '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/Homo_sapiens.GRCh38.93.sorted.gtf.gz'
print('loading gtf file')
gtf <- import(gtf_file, format = 'gtf')

utr3_gtf <- gtf[gtf$type=='three_prime_utr',]
length(unique(utr3_gtf$gene_id))
#19485
length(unique(utr3_gtf$transcript_id))
#74365

genes.Granges <- split(as.data.table(utr3_gtf),by='gene_id')
reduced.3utr <- do.call(rbind,lapply(genes.Granges, function(x){
  gene_id <- unique(x$gene_id)
  glist <- makeGRangesListFromDataFrame(x)
  reduced <- as.data.table(Reduce(union,glist))
  reduced$gene_id <- gene_id
  return(reduced)
}))

write.table(reduced.3utr,file='reduced.3utr.txt',sep='\t',col.names=T,row.names=F,quote=F)

GrList <- makeGRangesListFromDataFrame(reduced.3utr,
                                       keep.extra.columns = TRUE,
                                       split.field = 'gene_id')

all_width <- c()
for(i in 1:length(GrList)){
  all_width[i] <- sum(width(GrList[i]))
}
median(all_width)
#1585
df <- data.frame(length_3utr=all_width,gene_id=names(GrList))
write.table(df,file='utr_length.txt',col.names = T,row.names = F,quote = F,sep='\t')

pdf('../Analysis/figures/length_of_3utr_hs.pdf',8,6)
hist(log10(all_width),
     main= paste("Length of 3'UTR for Ensembl human genes",'median = 1585nt',
                 sep='\n'),
     ylab = 'Number of genes (total 19485)',
     xlab="Length of 3'UTR,log10(nt)",breaks = 20)
dev.off()




