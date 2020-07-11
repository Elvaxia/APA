setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP')
library(here)
library(data.table)
library(dplyr)
##k562
filenames <- list.files(here('K562'))
files <- list()
for ( i in filenames ) { files[[i]] <- fread(here('K562',i) , drop=1)}
names(files) <- gsub('IdentifiedPAs.','',names(files))
names(files) <- gsub('.txt','',names(files))
sample.info <- fread(here('RBP.sample.info.txt'),header=T)

####hepg2
hepg2.filenames <- list.files(here('ChangePoint/HepG2'))
hepg2.files <- list()
for ( i in hepg2.filenames ) { hepg2.files[[i]] <- fread(here('ChangePoint/HepG2',i) , drop=1)}
names(hepg2.files) <- gsub('IdentifiedPAs.','',names(hepg2.files))
names(hepg2.files) <- gsub('.txt','',names(hepg2.files))
sample.info.hepg2 <- fread(here('RBP.sample.HepG2.info.txt'))
###

##
ref.anno <- fread(here('src/reduced.3utr.txt'),header = T)
ref.anno.m <- ref.anno%>%group_by(gene_id)%>%mutate(st=min(start),ed=max(end))%>%
  dplyr::select(seqnames,strand,gene_id,st,ed)%>%unique()
###
filter.apa <- function(rbp.res){
  #rbp <- names(rbp.res)
  #rbp.res <- rbp.res[1][[1]]
  rbp <- gsub('_.*','',colnames(rbp.res)[1])
  #1.no APA
  rbp.res <- rbp.res[!is.na(PA),]
  #2.S_coverage lower than 10,convert RUD to be NA
  newfiles <- rbp.res
  newfiles$PA.ID <- paste(newfiles$gene_id,newfiles$PA,sep = '|')
  #remove PAS within 200bp distance to the start and end of region
  PAS <- merge(newfiles[,c('PA','gene_id','PA.ID')],ref.anno.m,by='gene_id',all.x = T)
  PAS <- PAS[abs(PAS$PA- PAS$st) >=200 & abs(PAS$PA-PAS$ed) >=200,]
  test.files <- merge(newfiles,PAS,by='PA.ID',all.y = T)
  test.files$PA.ID <- paste0(test.files$seqnames,':',test.files$PA.x,':',test.files$strand)
  #3.mean RUD of all samples greater than  1.
  #4.coverage of cUTR should be over 20
  rud.col.n <- grep('RUD',colnames(test.files))
  ruddt <- test.files[,..rud.col.n]
  mean_RUD <- rowMeans(ruddt,na.rm = T)
  sd_RUD <- apply(ruddt,1,function(x)sd(x,na.rm=T))
  median_RUD <- apply(ruddt,1,function(x)median(x,na.rm=T))
  N.NAs <- apply(ruddt,1,function(x)sum(is.na(x)))
  s.cov.col.n <- grep('S_Coverage',colnames(test.files))
  s.cov.dt <- test.files[,..s.cov.col.n]
  S.cov <- apply(s.cov.dt, 1, median)
  S.cov_na <- apply(s.cov.dt,1,function(x){
    sum(x<=10)
  })
  merged.genes <- apply(ruddt,1,function(x)sum(x>1.1,na.rm=T))
  filtered <- test.files[S.cov>=20&N.NAs <=2&
                         merged.genes <=2,]
  filtered$id <- paste(filtered$gene_id.x,
                      filtered$seqnames,
                      filtered$strand,
                      filtered$st,
                      filtered$PA.y,
                      filtered$ed,
                      sep = ':'
  )
  return(filtered)
}


#evaluate the statistical significance of RUD difference=======================
rbp.DE.apa <- function(filtered.res, test = 't.test'){
  filtered.res <- as.data.frame(filtered.res)
  kd.L.cov.mean <- filtered.res %>% dplyr::select(grep('_KD_[1-9].L_Coverage',colnames(filtered.res))) %>% rowMeans
  kd.L.cov.sd <- filtered.res %>% dplyr::select(grep('_KD_[1-9].L_Coverage',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
  kd.S.cov.mean <- filtered.res %>% dplyr::select(grep('_KD_[1-9].S_Coverage',colnames(filtered.res))) %>% rowMeans
  kd.S.cov.sd <- filtered.res %>% dplyr::select(grep('_KD_[1-9].S_Coverage',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
  kd.rud.mean <- filtered.res %>% dplyr::select(grep('_KD_[1-9].RUD',colnames(filtered.res))) %>% rowMeans
  kd.rud.sd <- filtered.res %>% dplyr::select(grep('_KD_[1-9].RUD',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
  ctrl.L.cov.mean <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].L_Coverage',colnames(filtered.res))) %>% rowMeans
  ctrl.L.cov.sd <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].L_Coverage',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
  ctrl.S.cov.mean <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].S_Coverage',colnames(filtered.res))) %>% rowMeans
  ctrl.S.cov.sd <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].S_Coverage',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
  ctrl.rud.mean <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].RUD',colnames(filtered.res))) %>% rowMeans
  ctrl.rud.sd <- filtered.res %>% dplyr::select(grep('_Ctrl_[1-9].RUD',colnames(filtered.res))) %>% apply(.,1,function(x)sd(x))
 
  if(test=='fisher.test'){
      fisher.test.mtx <- cbind(kd.L.cov.mean,kd.S.cov.mean,ctrl.L.cov.mean,ctrl.S.cov.mean)
    suppressWarnings(p_value <- apply(fisher.test.mtx,1,function(x){
      con.mtx <- matrix(x,nrow=2)
      p.v <- fisher.test(con.mtx)$p.value
      return(p.v)
    }))
  }
  if(test=='t.test'){
    t.test.mtx <- filtered.res%>%dplyr::select(grep('.RUD',colnames(filtered.res)))
    print(colnames(t.test.mtx)[1])
    t.test.mtx[t.test.mtx > 1] <- NA
    suppressWarnings(p_value <- apply(t.test.mtx,1,function(x){
      if(NA%in%x){
        p.v <- NA
      }else{
        p.v <- t.test(x[1:2],x[3:4])$p.value
      }     
      return(p.v)
    }))
}
  delta_RUD <- kd.rud.mean - ctrl.rud.mean
  log2FC <- log2(kd.rud.mean/ctrl.rud.mean)
  res <- as.data.table(cbind(p_value,delta_RUD,log2FC,kd.rud.mean,ctrl.rud.mean))
  res$id <- filtered.res$id
  #res <- res[p_value < 0.05& abs(delta_RUD)>0.1,]
  rbp <- gsub('_.*','',colnames(filtered.res)[2])
  #write.table(res,
  #            file = here('Analysis/rawdata',paste0(rbp,'.filtered.cov.RUD.txt')),
  #            col.names = T, row.names = F,quote = F)
  return(res)
}


###
all.apaevents.k562 <- lapply(files, filter.apa)
#some rbps have no control bam files
all.apaevents.k562 <- all.apaevents.k562[- which(as.numeric(lapply(all.apaevents.k562,ncol)) == 16)]

k562.de <- lapply(all.apaevents.k562,rbp.DE.apa)
new.k562.de <- k562.de

###
all.apaevents.hepg2 <- lapply(hepg2.files,filter.apa)
all.apaevents.hepg2 <- all.apaevents.hepg2[- which(as.numeric(lapply(all.apaevents.hepg2,ncol)) == 16)]
#some rbps have no control bam files
hepg2.de <- lapply(all.apaevents.hepg2,rbp.DE.apa)
##

length(hepg2.de)
load(here('Analysis/rawdata/rbp.testres.RData'))
save(new.k562.de,hepg2.de,file=here('Analysis/rawdata/rbp.testres.2.21.RData'))




#========================================================================
for (i in dt[[1]]$id){
  sashimi.plot(id=ENSG00000100614:14:+:60277088:60293189:60299087,rbp =rbp)
}
#===============











