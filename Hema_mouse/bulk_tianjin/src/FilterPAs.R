setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/')
source('./src/functions.hema_mm.R')

library(data.table)
library(dplyr)
filenames <- list.files('./ChangePoint/')
files <- list()
for(i in filenames){files[[i]] <- fread(paste0('./ChangePoint/',i),drop=1)}

files <- do.call(rbind,files)
#filtering out some events
#1.no APA
files <- files[!is.na(PA),]
#2.S_coverage lower than 10,convert RUD to be NA
sample.info <- fread('./Hema_sampleinfo.txt',header=T)
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
#3.mean RUD of all samples greater than  1.
#4.coverage of cUTR should be over 20

mean_RUD <- rowMeans(test.files[,70:103],na.rm = T)
sd_RUD <- apply(test.files[,70:103],1,function(x)sd(x,na.rm=T))
median_RUD <- apply(test.files[,70:103],1,function(x)median(x,na.rm=T))
N.NAs <- apply(test.files[,70:103],1,function(x)sum(is.na(x)))
S.cov <- apply(test.files[,36:69], 1, median)
S.cov_na <- apply(test.files[,36:69],1,function(x){
  sum(x<=10)
})
merged.genes <- apply(test.files[,70:103],1,function(x)sum(x>1.1,na.rm=T))


#以LT_HSC为标准，需要该基因在LT_HSC中表达，
S.cov.hsc <- apply(test.files[,70:71],1,function(x)sum(is.na(x)))
S.cov.mean.hsc<- apply(test.files[,70:71],1,function(x)mean(x,na.rm=T))
filtered <- test.files[S.cov>=20&N.NAs <=32&
                         merged.genes <=2,]
filtered$id <- paste(filtered$gene_id.x,
                     filtered$seqnames,
                     filtered$strand,
                     filtered$st,
                     filtered$PA.y,
                     filtered$ed,
                     sep = ':'
)

#write.table(filtered,file='./Analysis/rawdata/filtered.RUD.txt',col.names=T,row.names=F,quote=F)

#========================================================================
#sashimi.plot(id=filtered$id[69])


#weighted utr length=======================================================
library(data.table)
library(dplyr)
filtered <- fread('./Analysis/rawdata/filtered.RUD.txt')
UTR <- filtered[,c('gene_id.x','PA.y','seqnames','strand','st','ed','id')]
cUTR <- UTR %>% mutate(start = ifelse(strand=='+',st,PA.y ),
                       end = ifelse(strand=='+', PA.y,ed))
cUTR <- cUTR[,c(1,3,4,7,8,9)]
colnames(cUTR)[1] <- 'gene_id'

aUTR <- UTR %>% mutate(start = ifelse(strand=='+',PA.y,st ),
                       end = ifelse(strand=='+', ed,PA.y))
aUTR <- aUTR[,c(1,3,4,7,8,9)]
colnames(aUTR)[1] <- 'gene_id'

library(GenomicFeatures)
anno.UTR <- fread('src/reduced.3utr.txt')
setDT(aUTR)
setkey(aUTR,'gene_id','seqnames','strand','start','end')

setDT(cUTR)
setkey(cUTR,'gene_id','seqnames','strand','start','end')

setkey(anno.UTR,'gene_id','seqnames','strand','start','end')

cUTR.t <- foverlaps(cUTR,anno.UTR)
ce <- function(x){
  if(x[3]=='+'){
    st <- x[4]
    ed <- min(x[5],x[9])
  }else{
    st <- max(x[4],x[8])
    ed <- x[5]
  }
  return(c(st,ed))
}
test <- apply(cUTR.t, 1, ce)
cUTR.t <- cbind(cUTR.t,t(test))
cUTR.t <- cUTR.t[,-c(8,9)]
colnames(cUTR.t)[c(4:6,8:9)] <- c('utr_st','utr_ed','utr_width','start','end')
cUTR <- makeGRangesListFromDataFrame(cUTR.t,
                                     keep.extra.columns = TRUE,
                                     split.field = 'id')

aUTR.t <- foverlaps(aUTR,anno.UTR)
ces <- function(x){
  if(x[3]=='+'){
    st <- max(x[4],x[8])
    ed <- x[5]
  }else{
    st <- x[4]
    ed <- min(x[5],x[9])
  }
  return(c(st,ed))
}
test <- apply(aUTR.t, 1, ces)
aUTR.t <- cbind(aUTR.t,t(test))
aUTR.t <- aUTR.t[,-c(8,9)]
colnames(aUTR.t)[c(4:6,8:9)] <- c('utr_st','utr_ed','utr_width','start','end')

#write.table(aUTR.t,file = 'Analysis/rawdata/aUTR.txt',col.names = T,row.names = F,quote=F,sep='\t')
#write.table(cUTR.t,file = 'Analysis/rawdata/cUTR.txt',col.names = T,row.names = F,quote=F,sep='\t')

aUTR <- makeGRangesListFromDataFrame(aUTR.t,
                                     keep.extra.columns = TRUE,
                                     split.field = 'id')
width.aUTR <- lapply(aUTR, function(x)sum(width(x)))
#hist(log(as.numeric(width.aUTR)))
width.cUTR  <- lapply(cUTR, function(x)sum(width(x)))
#hist(log(as.numeric(width.cUTR)))



total.length <- as.numeric(width.cUTR) + as.numeric(width.aUTR)
#weighted_utr_length WUL
identical(names(width.cUTR),filtered$id)
filtered$utr.length<- as.numeric(total.length)

rud <- filtered[,70:101]
rud[rud>1] <- 0
rud.bulk <- list()
for(i in seq(1,32,2)){
  res <- rowMeans(rud[,i:(i+1)])
  rud.bulk[[i]] <- res
}

rud.bulk <- as.data.frame(do.call(cbind,rud.bulk))
colnames(rud.bulk) <- unique(gsub('_.*','',colnames(rud)))
colnames(rud) <- gsub('.RUD','',colnames(rud))

rud.bulk[rud.bulk >1] <- 1

WUL.L <-  rud.bulk *as.numeric(total.length)
WUL.S <- (1-rud.bulk) * as.numeric(width.cUTR)
WUL <- WUL.L + WUL.S
WUL[is.na(WUL)] <- 0
WUL[WUL > max(as.numeric(width.cUTR)) ] <- 0

#cor wul
ge <- scale(t(WUL))

#

#pearson correlation
cor.res <- cor(t(ge),method = 'pearson',use='pairwise.complete.obs')
library(pheatmap)
library(RColorBrewer)


my_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="rgb")(n=100)

breaks = generate_breaks(cor.res,length(colors),center = T)
pdf('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/Analysis/figures/cor.utr.length.pdf',8,8)
pheatmap(cor.res, main = 'Correlation of weighted utr length',
         color = my_palette)
dev.off()
#cor rud 
ge.rud <- scale(t(rud.bulk))
res.rud <- hcut(ge.rud, k = 4, stand = T)
cor.res.rud <- cor(t(ge.rud),method = 'pearson',use='pairwise.complete.obs')
#
pdf('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/Analysis/figures/cor.rucp.pdf',8,8)
pheatmap(cor.res.rud, main = 'Correlation of RUCP',
         color = my_palette)
dev.off()



#WUL and gene expression
WUL <- as.data.frame(WUL)
WUL$gene_id <- filtered$gene_id.x
WUL <- WUL[!duplicated(WUL$gene_id),]

wul.list <- list()
for(i in colnames(WUL)[2:16]){
  wul.list[[i]] <- WUL[,i] -WUL[,1]
}
wul.change <- as.data.frame(do.call(cbind,wul.list))
rownames(wul.change) <- WUL$gene_id

res.list <- readRDS('./Analysis/rawdata/DESeq.res.rds')
filtered.res <- list()
for(i in 1:15){
    te <- as.data.frame(res.list[[i]][,c(2,6)])
    te <- te[!is.na(te$padj)&te$padj <0.05,]
    te$gene <- rownames(te)
    #te <- te[,c(1,3)]
    filtered.res[[i]] <- te[,-2]
}
names(filtered.res) <- names(res.list)[1:15]
df.list <- list()
for(i in 1:15){
  use <- filtered.res[[i]]
  df <- merge(use,wul.change,by.x='gene',by.y='row.names',all=F)
  df <- df[,c(1:2,i+2)]
  colnames(df)[2:3] <- c('Gene','utr')
df <- df[df$Gene!='-Inf'&df$utr!='-Inf'&!is.na(df$utr)&!is.na(df$Gene),]
df <- df%>% mutate(change= ifelse(utr >=100,'Lengthened',ifelse(utr <= -100, 'Shortened','unchanged')))
df <- df%>% mutate(regulated = ifelse(Gene >0,'Up','Down'))
df.list[[i]]<- df
}

changed.genes <- list()
for(i in 1:15){
  changed.genes[[i]] <- as.character(df.list[[i]][df.list[[i]]$change!='unchanged',]$gene)
}
names(changed.genes) <- names(filtered.res)
#log10 utr
df.list <- lapply(df.list,function(x){
  tmp <- x%>%mutate(logutr = ifelse(utr>=0,log1p(utr),-log1p(abs(utr))))
  return(tmp)
})

names(df.list) <- names(filtered.res)
cor.test(df.list[[10]][df.list[[10]]$change!='unchanged',][,1],df.list[[10]][df.list[[10]]$change!='unchanged',][,2])

for(i in 1:15){
  ggplot(df.list[[i]],aes(x=Gene,y=logutr,fill=change)) + geom_point(aes(color=change),size=4)+
  labs(y="Difference of weighted 3'UTR length(log)", x= "log2 fold change of gene expression") +
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 22),
                  legend.position = 'right', legend.title = element_blank(),
                  legend.text =element_text(size=20),
                  axis.title.x = element_text(size=22),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=22),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
                  scale_y_continuous(labels=c(-9,-4.5,0,4.5,9),limits = c(-9, 9))
  ggsave(file=paste0('./Analysis/figures/GENEvsUTR/',names(df.list[i]),'gene.vs.utr.log.pdf'))
}

for(i in 1:15){
  ggplot(df.list[[i]],aes(x=Gene,y=utr,fill=change)) + geom_point(aes(color=change),size=4)+
  labs(y="Difference of weighted 3'UTR length", x= "log2 fold change of gene expression") +
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 22),
                  legend.position = 'right', legend.title = element_blank(),
                  axis.title.x = element_text(size=22),
                  legend.text =element_text(size=20),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=22),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
                  scale_y_continuous(labels=c(-4000,-2000,0,2000,4000),limits = c(-4000, 4000))
  ggsave(file=paste0('./Analysis/figures/GENEvsUTR/',names(df.list[i]),'gene.vs.utr.pdf'))
}



for(i in 1:15){
  df.list[[i]]%>%filter(change!='unchanged')%>%group_by(change)%>%dplyr::count(regulated)%>%as.data.frame()%>%
  ggplot(aes(x=change,y=n,fill=regulated)) + geom_bar(position = 'dodge', stat='identity')+
  labs(y='Number of genes')+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 22),
                  legend.position = 'top', legend.title = element_blank(),
                  axis.title.x = element_blank(),
                  legend.text =element_text(size=20),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=22),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#4393c3','#b2182b'))
  ggsave(file=paste0('./Analysis/figures/GENEvsUTR/',names(df.list[i]),'gene.vs.utr.up.down.pdf'))
}

#dotplot shows correlation between gene expression and delta RUCP================================

load('./Analysis/rawdata/DESeq.ntd.RData')
#LT-HSC vs ST-HSC
delta.rud <- function(x){
  LT_hsc <- filtered[,'LT-HSC_1.RUD']
  colunm <- paste0(x,'_1.RUD')
  tmp <- filtered[,..colunm]
  res <- as.data.frame(LT_hsc -tmp)
  res <- data.frame(gene=filtered$gene_id.x,delta_RUCP=res[,1])
  res <- res[!duplicated(res$gene),]
  res <- res[res$delta_RUCP!=0,]
  return(res)
}
delta.rucp.list <- list()
for(i in names(df.list)){
  delta.rucp.list[[i]] <- delta.rud(i)
}
delta.gene.list <- list()
for(i in names(df.list)){
  LT_hsc <- assay(ntd)[,'LT-HSC_1']
  colunm <- paste0(i,'_1')
  tmp <- assay(ntd)[,colunm]
  res <- LT_hsc -tmp
  res <- data.frame(gene=rownames(assay(ntd)),delta_geneExp=res)
  res <- res[res$delta_geneExp!=0,]
  delta.gene.list[[i]] <- res
}

#merge delta rucp and gene expression

rucp.gene <- list()
for(i in names(delta.gene.list)){
  rucp <- delta.rucp.list[[i]]
  gene <- delta.gene.list[[i]]
  res <- merge(rucp,gene,by='gene',all=F)
  rucp.gene[[i]] <- res
}

#correlation of wul and gene expression
rud <- filtered[,70:101]
rud[rud>1] <- 0
WUL.L <-  rud *as.numeric(total.length)
WUL.S <- (1-rud) * as.numeric(width.cUTR)
WUL <- WUL.L + WUL.S
WUL[is.na(WUL)] <- 0
WUL[WUL > max(as.numeric(width.cUTR)) ] <- 0

WUL$gene_id <- filtered$gene_id.x
WUL <- WUL[!duplicated(WUL$gene_id),]
WUL <- as.data.frame(WUL)
wul.gene <- list()
for(i in colnames(WUL)[1:32]){
  wul<- data.frame(gene=WUL$gene_id,WUL=WUL[,i])
  gene <- assay(dds)[,gsub('.RUD','',i)]
  gene <- data.frame(gene=rownames(assay(dds)),expression=gene)
  res <- merge(wul,gene,by='gene',all=F)
  res <- res[res$WUL!=0&res$expression>=10,]
  #res$WUL <- log10(res$WUL)
  wul.gene[[i]] <- res
}

#only use changed genes
#ST_HSC
rud.gene.sig.wul <- list()
for(i in names(changed.genes)){
  gene.df.col <- names(rud.gene)[grep(i ,names(rud.gene))][1]
  tmp <- rud.gene[[gene.df.col]]
  tmp <- tmp[tmp$gene%in%changed.genes[[i]],]
  rud.gene.sig.wul[[i]] <- tmp
}
for(i in names(rud.gene.sig.wul)){
  data <- rud.gene.sig.wul[[i]]
  setDT(data)
  colnames(data)[2:3]<- c("xc","yc")
  ggplot(data = data, aes(x = xc, y = yc)) +
    geom_point(size = .1)+
    stat_density2d(geom = "raster", aes(fill = ..density.., alpha = ..density..), contour = FALSE) +
    scale_fill_viridis_c(guide = FALSE) +
    scale_alpha_continuous(guide = "none", range = c(0, 1))+
    labs(x = "RUCP", y = "Gene expression",title = i,subtitle = paste("N =", nrow(data)))+
    #theme_classic() +
    geom_smooth(method = "lm", lwd = 0.5) +
stat_cor(data = data[, .(xc, yc)], method = "pearson",
             size = 6, col = "red",
             label.y.npc = "top", label.x.npc = "center") +
    theme(panel.background = element_blank(),
          plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 22),
          axis.line=element_line(),
          axis.text = element_text(size = 16), legend.position = "none")
ggsave(file=paste0(here('Analysis/figures/geneCORrud/'),i,'geneVSrud.sig.gene.pdf'))
}



tmp <- rud.gene[[1]]
tmp <- tmp[tmp$gene%in%changed.genes[[10]],]
cor.test(tmp$rud,tmp$expression)




library(ggpubr)



for(i in names(wul.gene)){
  data <- wul.gene[[i]]
  setDT(data)
  colnames(data)[2:3]<- c("xc","yc")
  ggplot(data = data, aes(x = xc, y = yc)) +
    geom_point(size = .1)+
    stat_density2d(geom = "raster", aes(fill = ..density.., alpha = ..density..), contour = FALSE) +
    scale_fill_viridis_c(guide = FALSE) +
    scale_alpha_continuous(guide = "none", range = c(0, 1))+
    labs(x = "Weighted UTR length", y = "Gene expression",title = i,subtitle = paste("N =", nrow(data)))+
    #theme_classic() +
    geom_smooth(method = "lm", lwd = 0.5) +
stat_cor(data = data[, .(xc, yc)], method = "pearson",
             size = 6, col = "red",
             label.y.npc = "top", label.x.npc = "center") +
    theme(panel.background = element_blank(),
          plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 22),
          axis.line=element_line(),
          axis.text = element_text(size = 16), legend.position = "none")
        ggsave(file=paste0(here('Analysis/figures/geneCORwul/'),i,'geneVSwul.pdf'))
}

##rud and gene expression
rud <- filtered[,70:101]
rud[rud>1] <- 0
rud$gene <- filtered$gene_id.x
rud <- rud[!duplicated(rud$gene),]
rud <- as.data.frame(rud)

rud.gene <- list()
for(i in colnames(rud)[1:32]){
  RUD<- data.frame(gene=rud$gene,rud=rud[,i])
  gene <- assay(ntd)[,gsub('.RUD','',i)]
  gene <- data.frame(gene=rownames(assay(ntd)),expression=gene)
  res <- merge(RUD,gene,by='gene',all=F)
  res <- res[res$rud!=0&res$expression>=2,]
  rud.gene[[i]] <- res
}

for(i in names(rud.gene)){
  data <- rud.gene[[i]]
  setDT(data)
  colnames(data)[2:3]<- c("xc","yc")
  ggplot(data = data, aes(x = xc, y = yc)) +
    geom_point(size = .1)+
    stat_density2d(geom = "raster", aes(fill = ..density.., alpha = ..density..), contour = FALSE) +
    scale_fill_viridis_c(guide = FALSE) +
    scale_alpha_continuous(guide = "none", range = c(0, 1))+
    labs(x = "RUCP", y = "Gene expression",title = i,subtitle = paste("N =", nrow(data)))+
    #theme_classic() +
    geom_smooth(method = "lm", lwd = 0.5) +
stat_cor(data = data[, .(xc, yc)], method = "pearson",
             size = 6, col = "red",
             label.y.npc = "top", label.x.npc = "center") +
    theme(panel.background = element_blank(),
          plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 22),
          axis.line=element_line(),
          axis.text = element_text(size = 16), legend.position = "none")
        ggsave(file=paste0(here('Analysis/figures/geneCORrud/'),i,'geneVSrud.pdf'))
}



####=====================================================================
#RUD swift across celltype,all celltype compare to LT-HSC
library(ggridges)
library(ggpubr)
mycol <- colorRampPalette(brewer.pal(12, "Set3"))(16)
rownames(rud.bulk) <- filtered$id
rud.bulk <- as.data.frame(rud.bulk)
RUD.change <- list()
for( i in colnames(rud.bulk)[2:16]){
  res <- log2(rud.bulk[,1] / (rud.bulk[,i])
  RUD.change[[i]] <- res
}
RUD.change <- as.data.frame(do.call(cbind, RUD.change))

RUD.change[is.na(RUD.change)] <- 0



pdf('./Analysis/figures/log2fd.pdf',width=8,height=6)
melt(RUD.change)%>%ggplot(aes(variable,value,fill=variable)) + 
                   geom_boxplot(outlier.colour=NA) + 
                   labs(y = "Log2 of RUCP ratio :LT-HSC vs other celltypes") +
                   theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 16),
                  legend.position = 'none', legend.title = element_blank(),
                  axis.text.x = element_text(size=14,angle= -30),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=15),
                  axis.text.y  = element_text(size=15),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) +
                    scale_fill_manual(values=mycol)+
                    coord_cartesian(ylim = c(-2,2)) 
                    
                    #stat_compare_means(method = "anova", label.y = 2) + # Add global p-value
                    #stat_compare_means(aes(label = ..p.signif..),method = "t.test",
                    #label.y = 1.8, ref.group = "ST-HSC")
dev.off()

##shortend and lengthend
length.change.ery <- list()
for( i in colnames(rud.bulk)[15:16]){
  res <- as.numeric(rud.bulk[,14] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.ery[[i]] <- table(res)
}
length.change.ery <- do.call(cbind, length.change.ery)
length.change.lym <- list()
for( i in colnames(rud.bulk)[9:13]){
  res <- as.numeric(rud.bulk[,9] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.lym[[i]] <- table(res)
} 

length.change.lym <- do.call(cbind, length.change.lym)

length.change.all <- cbind(length.change.lym,length.change.ery)
colnames(length.change.all) <- c('B-CLP','CD4-CLP','CD8-CLP','NK-CLP','EryA-MEP','EryB-MEP')

pdf('./Analysis/figures/lengthchange.pdf',width=8,height=6)
melt(length.change.all)%>%ggplot(aes(Var2,value,fill=Var1)) + 
                   geom_bar(position = 'dodge', stat='identity') +
                   labs(y = "Nuber of genes shortend/lengthed") +
                   theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 16),
                  legend.position = 'top', legend.title = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x =element_text(size=15),
                  axis.title.y = element_text(size=15),
                  axis.text.y  = element_text(size=15),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) 
                    #stat_compare_means(method = "anova", label.y = 2) + # Add global p-value
                    #stat_compare_means(aes(label = ..p.signif..),method = "t.test",
                    #label.y = 1.8, ref.group = "ST-HSC")
dev.off()


####lengthchange all==========
length.change.ery <- list()
for( i in colnames(rud.bulk)[15:16]){
  res <- as.numeric(rud.bulk[,14] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.ery[[i]] <- table(res)
}
length.change.ery <- do.call(cbind, length.change.ery)

length.change.lym <- list()
for( i in colnames(rud.bulk)[10:13]){
  res <- as.numeric(rud.bulk[,9] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.lym[[i]] <- table(res)
} 
length.change.lym <- do.call(cbind, length.change.lym)


length.change.mye <- list()
for( i in colnames(rud.bulk)[6:8]){
  res <- as.numeric(rud.bulk[,5] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.mye[[i]] <- table(res)
} 
length.change.mye <- do.call(cbind, length.change.mye)


length.change.pro <- list()
for( i in colnames(rud.bulk)[c(2:5,9)]){
  res <- as.numeric(rud.bulk[,5] - rud.bulk[,i])
  res <- res[res!=0]
  res <- ifelse(res >0,'Shortened','Lengthend')
  length.change.pro[[i]] <- table(res)
} 
length.change.pro <- do.call(cbind, length.change.pro)




length.change.all <- cbind(length.change.lym,length.change.ery,length.change.mye,length.change.pro)
colnames(length.change.all) <- c('B-CLP','CD4-CLP','CD8-CLP','NK-CLP','EryA-MEP','EryB-MEP',
'Gn-GMP','Mo-GMP','Mp-GMP','ST_HSC-LT_HSC','MPP-LT_HSC','CMP-LT_HSC','CLP-LT_HSC')
rownames(length.change.all) <- c('Lengthened','Shortened')

pdf('./Analysis/figures/lengthchange.pdf',width=8,height=6)
melt(length.change.all)%>%ggplot(aes(Var2,value,fill=Var1)) + 
                   geom_bar(position = 'dodge', stat='identity') +
                   labs(y = "Nuber of genes shortened/lengthened") +
                   theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 16),
                  legend.position = 'top', legend.title = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x =element_text(size=15,angle= 90),
                  axis.title.y = element_text(size=15),
                  axis.text.y  = element_text(size=15),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) 
                    #stat_compare_means(method = "anova", label.y = 2) + # Add global p-value
                    #stat_compare_means(aes(label = ..p.signif..),method = "t.test",
                    #label.y = 1.8, ref.group = "ST-HSC")
dev.off()







