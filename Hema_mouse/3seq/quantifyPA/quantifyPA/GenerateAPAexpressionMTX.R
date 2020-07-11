#step1:group PA sites
#step2:load apa expression and filter express

setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/quantifyPA/all.pa')
library(data.table)
all.pa <- fread('all.PA.combined.txt')
colnames(all.pa) <- c('id','counts')
all.pa <- all.pa[,lapply(.SD, sum),by=id]

all.pa$chr <- do.call(c,lapply(strsplit(all.pa$id,'\\:'),function(x)x[1]))
all.pa$strand <- do.call(c,lapply(strsplit(all.pa$id,'\\:'),function(x)x[3]))
all.pa$pos <- as.numeric(do.call(c,lapply(strsplit(all.pa$id,'\\:'),function(x)x[2])))
setkey(all.pa,chr,strand,pos)

#step1:group PA

#difference between sites
all.PA <- all.pa

all.PA$dif<-c(0,diff(all.PA[,pos]))

ts <- all.PA
n.id<-ts[1,pos]
for (i in 2:nrow(ts)) {
  if (ts[i,dif]>=0 & ts[i,dif]<=24) {
    for (j in 1:i) {
      print(c(i,j))
      if (ts[i-j,dif]>=0 & ts[i-j,dif]<=24) {
        print("continue")
      } else {
        add.id<-ts[i-j,pos]
        break
      }
    }
  } else {
    add.id<-ts[i,pos]
  }
  n.id<-c(n.id,add.id)
}

all.PA$new.id<-n.id


PA.group <- data.frame(PA.ori= paste(all.PA$chr,all.PA$pos,
                                     all.PA$strand,sep=':'),
                       PA.grouped=paste(all.PA$chr,all.PA$new.id,
                                        all.PA$strand,sep=':'))

write.table(PA.group,file='PA.grouped.txt',sep='\t',
            col.names = T,row.names = F,quote = F)




#step2:
#load APA expression
PA.group <- read.table('PA.grouped.txt',header=T)

mypath <- '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/quantifyPA/all.pa/groupPA'
multmerge = function(mypath,id.group){
  #need to change the pattern accordingly !!!
  filenames = list.files(path=mypath, pattern=".txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    #for x add group id
    print(x)
    ori <- fread(x)
    colnames(ori) <- c('PA.ori','reads')
    tmp<-merge(ori,id.group,by="PA.ori")
    setkey(tmp,PA.grouped)
    tmps<-tmp[,lapply(.SD, sum),by=PA.grouped, .SDcols = 'reads']
    setkey(tmps,PA.grouped)
    colnames(tmps)[2] <- gsub('/','',gsub('.PA.txt','',
                                          gsub(mypath,'',x)))
    return(tmps)})
  Reduce(function(x,y) {merge(x,y,all=T,by="PA.grouped")}, datalist)
}

Expr <- multmerge(mypath = mypath, id.group = PA.group)
rs <- rowSums(Expr[,-1],na.rm=T)
nas <- apply(Expr[,-1],1,function(x)sum(is.na(x)))
Expr1 <- Expr[rs>10&nas <51,]

PAs <- as.character(Expr1$PA.grouped)
write.table(PAs,file='PA.grouped.filtered.txt',col.names = F,row.names = F,quote = F)
saveRDS(Expr,file='PA.Expr.ori.rds')
saveRDS(Expr1,file='PA.Expr.filter.rds')


