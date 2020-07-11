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

#opt$gtf <- '/mnt/raid61/data/Hematopoiesis/Mouse/Mus_musculus.GRCm38.95.gtf.gz'
gtf_file <- opt$gtf


gtf_file <- '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz'
print('loading gtf file')
gtf <- import(gtf_file, format = 'gtf')

utr3_gtf <- gtf[gtf$type=='three_prime_utr',]
length(unique(utr3_gtf$gene_id))
#21090
length(unique(utr3_gtf$transcript_id))
#50431

genes.Granges <- split(as.data.table(utr3_gtf),by='gene_id')
reduced.3utr <- do.call(rbind,lapply(genes.Granges, function(x){
  gene_id <- unique(x$gene_id)
  glist <- makeGRangesListFromDataFrame(x)
  reduced <- as.data.table(Reduce(union,glist))
  reduced$gene_id <- gene_id
  return(reduced)
}))

write.table(reduced.3utr,file='reduced.3utr.txt',sep='\t',col.names=T,row.names=F)

GrList <- makeGRangesListFromDataFrame(reduced.3utr,
                                       keep.extra.columns = TRUE,
                                       split.field = 'gene_id')

all_width <- c()
for(i in 1:length(GrList)){
  all_width[i] <- sum(width(GrList[i]))
}
df <- data.frame(length_3utr=all_width,gene_id=names(GrList))
write.table(df,file='utr_length.txt',col.names = T,row.names = F,quote = F,sep='\t')

pdf('../Analysis/figures/length_of_3utr_mouse.pdf',8,6)
hist(log10(all_width),
     main= paste("Length of 3'UTR for Ensembl mouse genes",'median = 1164nt',
                 sep='\n'),
     ylab = 'Number of genes (total 21090)',
     xlab="Length of 3'UTR,log10(nt)",breaks = 20)
dev.off()











txdb <- makeTxDbFromGRanges(gtf, drop.stop.codons=FALSE)
exons <- exonsBy(txdb, use.name=TRUE)
introns <- intronsByTranscript(txdb, use.name=TRUE)
utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
cds <- cdsBy(txdb, "tx", use.names=TRUE)

##Extract 3UTR ---------------------------------
# Selcet 3UTR on the terminal exon,

# Remove the 5UTR on the same exon
# Remove artifacts: antisense RNA
# ===NEED TODO:after compute the UTR3 usage,=====
# ===filter out the expressed overlapped genes==

Extend3utr <- function(utr3, exons, utr5){
  utr3.exon <- as.data.table(utr3) %>% group_by(group_name) %>%
    filter(exon_rank == max(exon_rank)) 
  last.exon <- as.data.table(exons) %>% group_by(group_name) %>%
    filter(exon_rank == max(exon_rank))
  utr3.extracted <- merge(utr3.exon,last.exon,
                          by = 'group_name', all.y =FALSE) %>%
    mutate(start = ifelse(strand.x == '+',
                          start.y, start.x),
           end = ifelse(strand.x == '+',
                        end.x,end.y))
  
  utr3.extracted <- utr3.extracted[,-c(2,8,9,11,12,16,17,19)]
  colnames(utr3.extracted) <- c('enst.id', 'chr', 'utr3.start', 'utr3.end', 
                                'utr3.width', 'strand', 'exon.rank',
                                'exon.start', 'exon.end', 'exon.width', 'exon_name',
                                'start','end')
  #resize exon with 5'utr
  utr5.dt <- as.data.table(utr5)
  remove <- merge(utr3.extracted, utr5.dt, by = 'exon_name',
                  all.x = TRUE, all.y = FALSE) %>%
    filter(is.na(seqnames) | seqnames == chr) %>% #chrX has same exon_name as chrY
    mutate(start.x = ifelse(strand.x == '+', 
                            ifelse(is.na(end.y), start.x, end.y + 1), start.x),
           end.x = ifelse(strand.x == '-', 
                          ifelse(is.na(start.y), end.x, start.y - 1), end.x))
  remove <- remove[,1:13]
  colnames(remove)[c(7,12,13)] <- c('strand','start','end')
  return(remove)
}

AnnotateUTR3 <- function(utr3, exons, utr5, gtf){
  utr3.ref <- Extend3utr(utr3, exons, utr5)
  
  gene.id <- as.data.table(gtf[, c('gene_id', 'gene_biotype', 'gene_name',
                                   'transcript_id', 'transcript_biotype')])[, c(2,3,6:10)] %>%
    .[!duplicated(transcript_id), ] %>% filter(!is.na(transcript_id))
  
  utr.ref <- merge(utr3.ref, gene.id, by.x = 'enst.id', by.y = 'transcript_id',
                   all.x = TRUE, all.y = FALSE) %>%
    filter(gene_biotype != 'antisense_RNA')
  #remove artifacts:antisense RNA
  colnames(utr.ref)[12:15] <- c('start', 'end', 'tx.start', 'tx.end')
  print(paste(length(unique(utr.ref$gene_id)), 'genes',
              length(unique(utr.ref$enst.id)), 'transcripts'))
  return(utr.ref)
  
}


print('extract 3utr')

utr3.extracted <- AnnotateUTR3(utr3, exons, utr5, gtf)

library(SummarizedExperiment)
#TODO: After compute the usage,remove the expressed overlapped genes
write.table(utr3.extracted, file = paste(opt$outdir, 'UTR3.ref.txt',sep='/'),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')


#write.table(utr3.extracted, file = 'UTR3.ref.txt',
#            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')



##Extract intron-----------------------------------------
# Select introns in CDS.
# Discard introns range under 2000bp.
# Focuse on IpA of protein coding genes,
# remove introns that potentially originated from microRNAs, 
# small nucleolar RNAs, and retrotransposons.



ExtractIntron <- function(introns, cds, gtf){
  cds.combine <- as.data.table(cds) %>%
    group_by(group_name) %>%
    mutate(start = min(start), end = max(end)) %>%
    mutate(width = end -start) %>%
    filter(!duplicated(group_name)) %>% 
    setDT(key = c('seqnames', 'strand', 'start', 'end'))
  
  intron <- as.data.table(introns)[,-1] %>% filter(width >= 2 * 1e3) %>%
    setDT(key = c('seqnames', 'strand', 'start', 'end'))
  
  intron.cds <- foverlaps(cds.combine, intron, nomatch = 0) %>%
    filter(group_name == i.group_name) %>% .[,c(1:6,9:11)]
  
  colnames(intron.cds) <- c('chr', 'strand', 'transcript_id',
                            'intron.start', 'intron.end', 'intron_width',
                            'tx.cds.start', 'tx.cds.end', 'tx.cds.width')
  
  gene.id <- as.data.table(gtf[, c('gene_id', 'gene_biotype', 'gene_name',
                                   'transcript_id', 'transcript_biotype')])[, c(2,3,6:10)] %>%
    .[!duplicated(transcript_id), ] %>% filter(!is.na(transcript_id))
  
  intron.ref <- merge(intron.cds, gene.id, by = 'transcript_id',
                      all.x = TRUE, all.y = FALSE) %>%
    filter(transcript_biotype == 'protein_coding') 
  colnames(intron.ref)[c(1,4,5,10:11)] <- c('enst.id', 'start', 'end','tx.start', 'tx.end')
  return(intron.ref)
}

print('extracting intron')
intron.extracted <- ExtractIntron(introns, cds,  gtf)

write.table(intron.extracted, file = paste(opt$outdir, 'INTRON.ref.txt',sep='/'),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')


#write.table(intron.extracted, file = 'INTRON.ref.txt',
 #           quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

##------------------------------------------------
## criteria of filter and reassign the events 

#Some genes in the genome overlap with each
#other. such genes  were removed from further analysis.

#we were interested in investigating the IpA isoforms of protein coding genes, 
#events falling in introns that potentially originated from microRNAs, 
#small nucleolar RNAs, and retrotransposons were also removed

#peaks in introns that were close to the end of an 
#opposite strand 3ʹUTR were also removed

#There are genes where the end of the 3ʹUTR might fall
#in the intron of the downstream gene on the same strand.
#This would also create peaks in introns that are contributed by the preceding gene
#the intron  within 5000 nt of the 3ʹ end of the 3ʹUTR were reassigned to 3UTR

##---------------------------------------------
## compute usage
##the usage of IpA isoforms was calculated 

#we filtered for
#robustly expressed isoforms by imposing TPM and usage cutoffs
