#annotate PA site
#input PA site format:chr:site:strand
#output one file
suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(annotatr)
  library(dplyr)
})

site <-'../quantifyPA/all.pa/PA.grouped.filtered.txt'
txt2bed <-function(file){
  df <- fread(file,header = F)
  df <- as.data.frame(do.call(rbind, strsplit(as.character(unique(df$V1)),":")))
  colnames(df) <- c("chr","end","strand")
  df$chr <- gsub('chr','',df$chr)
  df$start <- as.integer(as.character(df$end)) - 1
  df$score <- '.'
  df <- makeGRangesFromDataFrame(df)
  return(df)
}
site.bed <- txt2bed(site)

gtf_file <- '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz'
gtf <- import(gtf_file, format = 'gtf')

#we consider PA sites withn 24nt to annotated site,so extended the 3'end of gtf
gtf.new <- as.data.frame(gtf) %>% 
  mutate(start=ifelse(strand=='+',start,start-24),
                                         end=ifelse(strand=='+',end+24,end))
gtf.new <- makeGRangesFromDataFrame(gtf.new,keep.extra.columns=T)
hits.res <- findOverlaps(site.bed,gtf.new)
site.bed.queryed <- as.data.table(site.bed[hits.res@from,])
gtf.queryed <- as.data.table(gtf.new[hits.res@to,])
site.res <- cbind(site.bed.queryed,gtf.queryed)
colnames(site.res)[1:5] <- c('chr','st','PAs','sitewidth','strand.ps')
site.res$PA <- paste(site.res$chr,site.res$PAs,site.res$strand.ps,sep=':')

UTR3 <- site.res[type=='three_prime_utr',]
UTR3.pa <- unique(UTR3$PA)
UTR5 <- site.res[type=='five_prime_utr',]
UTR5.pa <- unique(UTR5$PA)
CDS <- site.res[type=='CDS',]
CDS.pa <- unique(CDS$PA)

##extend 2kb of gtf to find extended 3utr
extended2kb.3utr <- as.data.frame(gtf) %>% 
  mutate(start=ifelse(strand=='+',start,start-2000),
         end=ifelse(strand=='+',end+2000,end))
extended2kb.3utr <- makeGRangesFromDataFrame(extended2kb.3utr,keep.extra.columns=T)
  
hits.res <- findOverlaps(site.bed,extended2kb.3utr)
site.bed.queryed <- as.data.table(site.bed[hits.res@from,])
gtf.queryed <- as.data.table(extended2kb.3utr[hits.res@to,])
site.res <- cbind(site.bed.queryed,gtf.queryed)
colnames(site.res)[1:5] <- c('chr','st','PAs','sitewidth','strand.ps')
site.res$PA <- paste(site.res$chr,site.res$PAs,site.res$strand.ps,sep=':')

UTR3 <- site.res[type=='three_prime_utr',]
UTR3.2kb.pa <- unique(UTR3$PA)

genes <- site.res[type=='gene',]
genes.pa <- unique(genes$PA)

#annotated PA type priority:
#3utr,3utr_2kb,cds,5utr,intron,
all.site <- fread(site,header = F)$V1
final.3utr.pa <- UTR3.pa
final.3utr.2kb.pa <- UTR3.2kb.pa[!UTR3.2kb.pa%in%UTR3.pa]
final.CDS <- CDS.pa[!CDS.pa%in%c(final.3utr.pa,final.3utr.2kb.pa)]
final.5utr <- UTR5.pa[!UTR5.pa%in%c(final.3utr.pa,final.3utr.2kb.pa,final.CDS)]
final.intron <- genes.pa[!genes.pa%in%c(final.3utr.pa,final.3utr.2kb.pa,
                                      final.CDS,final.5utr)]

final.intergenic <- all.site[!all.site%in%c(final.3utr.pa,final.3utr.2kb.pa,
                                           final.CDS,final.5utr,final.intron)]


##combine results
final.pa <- site.res[site.res$type=='gene',]
final.pa$PA.type <- 'INTERGENIC'
final.pa[PA%in%final.3utr.pa]$PA.type <- '3_UTR'
final.pa[PA%in%final.3utr.2kb.pa]$PA.type <- '3_UTR_extended'
final.pa[PA%in%final.CDS]$PA.type <- 'CDS'
final.pa[PA%in%final.5utr]$PA.type <- '5_UTR'
final.pa[PA%in%final.intron]$PA.type <- 'INTRON'
final.pa <- final.pa[,-c(4:9,11:14,16:18,22:32)]
colnames(final.pa)[2:3] <- c('start','end')

interpa <- as.data.frame(do.call(rbind, strsplit(as.character(final.intergenic),":")))
colnames(interpa) <- c("chr","end","strand")
interpa$start <- as.integer(as.character(interpa$end)) - 1
interpa$PA <- final.intergenic
interpa$PA.type <- 'INTERGENIC'


final.pa <- rbind(final.pa,interpa, fill=T)


write.table(final.pa,file='annotated.PA.txt',
            col.names = T,row.names = F,sep='\t',quote=F)




