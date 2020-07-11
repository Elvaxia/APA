
library(data.table)
library(dplyr)
filenames <- list.files('./Changepoint/')
files <- list()
for(i in filenames){files[[i]] <- fread(paste0('./Changepoint/',i),drop=1)}

files <- do.call(rbind,files)
#filtering out some events
#1.no APA
files <- files[!is.na(PA),]
#2.S_coverage lower than 10,convert RUD to be NA
sample.info <- fread('hs_hema_sampleinfo.txt',header=T)
samples <- sample.info$sample
newfiles <- files
newfiles$PA.ID <- paste(newfiles$gene_id,newfiles$PA,sep = '|')


ref.anno <- fread('./src/reduced.3utr.txt',header = T)
#remove PAS within 150bp distance to the start and end of region
ref.anno.m <- ref.anno%>%group_by(gene_id)%>%mutate(st=min(start),ed=max(end))%>%
  select(seqnames,strand,gene_id,st,ed)%>%unique()

PAS <- merge(newfiles[,c('PA','gene_id','PA.ID')],ref.anno.m,by='gene_id',all.x = T)

dist_s <- abs(PAS$PA- PAS$st)
dist_e <- abs(PAS$PA-PAS$ed)

PAS <- PAS[dist_s >=200 &dist_e >=200,]

test.files <- merge(newfiles,PAS,by='PA.ID',all.y = T)
test.files$PA.ID <- paste0(test.files$seqnames,':',test.files$PA.x,':',test.files$strand)
#3.mean RUD of all samples greater than  1

mean_RUD <- rowMeans(test.files[,130:193],na.rm = T)
sd_RUD <- apply(test.files[,130:193],1,function(x)sd(x,na.rm=T))
median_RUD <- apply(test.files[,130:193],1,function(x)median(x,na.rm=T))
N.NAs <- apply(test.files[,130:193],1,function(x)sum(is.na(x)))


filtered <- test.files[median_RUD < 1 &N.NAs <=30,]
filtered$PA.ID <- paste0(filtered$seqnames,':',filtered$PA.x,':',filtered$strand)
length(unique(filtered$gene_id.x))
#7013
write.table(filtered,file='./Analysis/rawdata//filtered.RUD.txt',col.names=T,row.names=F,quote=F)




