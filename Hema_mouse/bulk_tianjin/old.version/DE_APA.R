#setwd('./all_chrs')
library(data.table)
chrs <- c(seq(1,19),'X','Y')
files <- list()
for(i in chrs){files[[i]] <- fread(paste0(i,'/IdentifiedPAs.all.txt'),header=T)}
files <- do.call(rbind,files)
#filtering out some events
#1.no APA
files <- files[!is.na(PA),]
#2.S_coverage lower than 10,convert RUD to be NA
sample.info <- fread('../Hema_sampleinfo.txt',header=T)
samples <- sample.info$sample
newfiles <- files
for(s in samples){
  for(i in 1:nrow(files)){
    print(c(s,i))
    s_c <- paste0(s,'.S_Coverage')
    rud <- paste0(s,'.RUD')
    if(files[i,..s_c] < 10){
      newfiles[i,paste0(s,'.RUD')] <- NA
    }
  }
}
write.table(newfiles, file='../DE_APA/IdentifiedPAs.all.chrs.cutoffCoverage.txt',col.names=T,
            row.names=F,sep='\t',quote=F)

newfiles <- files
for(s in samples){
  for(i in 1:nrow(files)){
    print(c(s,i))
    s_c <- paste0(s,'.S_Coverage')
    rud <- paste0(s,'.RUD')
    if(files[i,..s_c] < 50){
      newfiles[i,paste0(s,'.RUD')] <- NA
    }
  }
}
write.table(newfiles, file='../DE_APA/IdentifiedPAs.all.chrs.cutoffCoverage50.txt',col.names=T,
            row.names=F,sep='\t',quote=F)






newfiles <- fread('./DE_APA/IdentifiedPAs.all.chrs.cutoffCoverage50.txt',header = T)
newfiles$PA.ID <- paste(newfiles$ID,newfiles$PA,sep = '|')


ref.anno <- fread('GRCm38.3utr.unextended.bed',header = F)
#remove PAS within 150bp distance to the start and end of region
colnames(ref.anno) <- c('chr','start','end','ID','score','strand')

PAS <- merge(newfiles[,c('PA','ID','PA.ID')],ref.anno,by='ID',all.x = T)

dist_s <- abs(PAS$PA- PAS$start)
dist_e <- abs(PAS$PA-PAS$end)

PAS <- PAS[dist_s >=200 &dist_e >=200,]

test.files <- merge(PAS,newfiles,by='PA.ID',all.x=T)
#3.mean RUD of all samples greater than  1

mean_RUD <- rowMeans(test.files[,77:110],na.rm = T)
sd_RUD <- apply(test.files[,77:110],1,function(x)sd(x,na.rm=T))
median_RUD <- apply(test.files[,77:110],1,function(x)median(x,na.rm=T))
N.NAs <- apply(test.files[,77:110],1,function(x)sum(is.na(x)))


filtered <- test.files[median_RUD < 1 &N.NAs <=30,]
write.table(filtered,file='./DE_APA/filtered.RUD.txt',col.names=T,row.names=F,quote=F)

##PA sites detected in 3'Seq

PA.3seq <- fread('/mnt/raid61/APA/Project/HSC_mm/mm_apa/basedistr/allPA.txt',
                 header = F)
PA.3seq$chr <- do.call(c,lapply(strsplit(PA.3seq$V1,'\\:'),function(x)x[1]))
PA.3seq$start <- as.numeric(do.call(c,lapply(strsplit(PA.3seq$V1,'\\:'),
                                             function(x)x[2]))) -24
PA.3seq$end <- PA.3seq$start +48
PA.3seq$strand <- do.call(c,lapply(strsplit(PA.3seq$V1,'\\:'),function(x)x[3]))

PA.3seq$seq3.id <- PA.3seq$V1
annoDT <- fread('/mnt/raid61/APA/Project/mm10.annotated.PA.database.txt',
                header=T)
annoDT$chr <- gsub('chr','',annoDT$chr)
annoDT$anno.id <- paste(annoDT$chr,annoDT$start,annoDT$strand,sep = ':')

setkey(annoDT,strand,chr,start,end)
setkey(PA.3seq,strand,chr,start,end)



PA.bulk <- data.table(PA = filtered$PA.y,ID=filtered$ID.y)
PA.bulk$start <- PA.bulk$PA - 100
PA.bulk$chr <- do.call(c,lapply(strsplit(PA.bulk$ID,'\\|'),function(x)x[3]))
PA.bulk$end <- PA.bulk$PA +100
PA.bulk $strand <- do.call(c,lapply(strsplit(PA.bulk$ID,'\\|'),function(x)x[4]))
PA.bulk$bulk.id <- paste0(PA.bulk$chr,':',PA.bulk$PA,':',PA.bulk$strand)
setkey(PA.bulk,strand,chr,start,end)
#PA.3seq and annoDT
over<-foverlaps(PA.3seq,annoDT,nomatch=0)
over1<- foverlaps(PA.3seq,PA.bulk,nomatch = 0)
over2 <- foverlaps(annoDT,PA.bulk,nomatch=0)


use <- unique(c(unique(over1$bulk.id),unique(over2$bulk.id)))
length(unique(over1$bulk.id))
length(unique(over2$bulk.id))
length(use)/length(unique(PA.bulk$bulk.id))
#
length(unique(over2$bulk.id))/length(unique(PA.bulk$bulk.id))
use1 <- unique(over1$bulk.id)[!unique(over1$bulk.id)%in%unique(over2$bulk.id)]
PA.bulk <- PA.bulk%>%mutate(new.PA = ifelse(strand=='+',PA,PA+50))
PA.bulk$new.bulk.id <- paste0(PA.bulk$chr,':',PA.bulk$new.PA,':',PA.bulk$strand)

write.table(unique(PA.bulk$new.bulk.id),
            file='/mnt/raid61/APA/Project/Hema_mouse/MEME/PA.txt',
            col.names=F,row.names=F,quote = F)

anno_and_3seq <- intersect(unique(over1$bulk.id),unique(over2$bulk.id))
only_3seq <- fsetdiff(data.table(unique(over1$bulk.id)),data.table(unique(over2$bulk.id)))$V1
only_anno <- fsetdiff(data.table(unique(over2$bulk.id)),data.table(unique(over1$bulk.id)))$V1

x<-c(length(anno_and_3seq),length(only_3seq),length(only_anno))

labels <-  c("Annotated by database and 3'seq ","only annotated by 3'seq","only annotated by database","This Study")
piepercent<- paste0(paste0(round(100*x/sum(x), 1),"%"),"(",x,")")
plot<-pie(x, labels = piepercent, 
          main = paste("PA sites detected by Microwell-Seq within",range,"bp",sep
                                               =""),cex=1.3)
legend("topleft", c("only Gencode","Gencode and PolyA_DB","only PolyA_DB","This Study"), cex = 1,
       fill = color, bty = "n")




