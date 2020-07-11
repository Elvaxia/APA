#DE_APA
library(data.table)
files <- fread('IdentifiedPAs.all.txt',header=T)
sample.info <- fread('../ALL_sample.info.txt',header=T)
samples <- sample.info$sample
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
write.table(newfiles, file='../IdentifiedPAs.all.chrs.cutoffCoverage.txt',col.names=T,
            row.names=F,sep='\t',quote=F)


ref <- fread('../ref/GRCh38.3utr.unextended.bed')
colnames(ref) <- c('chr','start','end','ID','score','strand')
sample.info <- fread('../ALL_sample.info.txt',header=T)
all.res <- fread('../IdentifiedPAs.all.chrs.cutoffCoverage.txt',header=T)
N.NAs <- apply(all.res[,41:60],1,function(x)sum(is.na(x)))
median.psi <- apply(all.res[,41:60],1,function(x)median(x,na.rm = T))
filtered <- all.res[N.NAs <= 18,,41:60]

files <- filtered

samples <- seq(1,10,1)
samples.mtx <- list()
for(i in samples){
  col <- c(paste0('ALL-',i,'.L_Coverage'),
           paste0('ALL-',i,'.S_Coverage'),
           paste0('CTRL-',i,'.L_Coverage'),
           paste0('CTRL-',i,'.S_Coverage'),
           paste0('ALL-',i,'.RUD'),
           paste0('CTRL-',1,'.RUD'))
  dt <- files[,..col]
  suppressWarnings(p_value <- apply(dt,1,function(x){
    con.mtx <- matrix(as.numeric(x[1:4]),nrow=2)
    p.v <- fisher.test(con.mtx)$p.value
    return(p.v)
  }))
  FDR <- p.adjust(p_value)
  delta_RUD <- dt[,5]-dt[,6]
  log2FC <- log2(dt[,5]/dt[,6])
  res <- cbind(p_value,FDR,delta_RUD,log2FC)
  colnames(res) <- c('p_value','FDR','delta_RUD','log2FC')
 samples.mtx[[i]] <- res
}

saveRDS(samples.mtx,file='samples.mtx.rds')

samples.mtx <- readRDS('samples.mtx.rds')

test.res <- function(x){
  x$ID <- paste(files$ID,files$PA,sep='|')
  filtered.res <- x[FDR <=0.05& abs(delta_RUD) >= 0.2&abs(log2FC)>=0.5,]$ID
  return(filtered.res)
}


all.filtered.id <- do.call(c,lapply(samples.mtx, test.res))
genes <- do.call(c,lapply(strsplit(all.filtered.id,'\\|'), function(x)x[1]))

all.filtered.id <- all.filtered.id[table(all.filtered.id) >=2]

all.filtered.id <- unique(all.filtered.id)

files$PA.ID <- paste(files$ID,files$PA,sep='|')
test.all.pa <- files[PA.ID%in%unique(all.filtered.id),]
test.all.pa$mean_RUD_ALL <- rowMeans(test.all.pa[,41:50],na.rm = T)
test.all.pa$mean_RUD_CTRL <- rowMeans(test.all.pa[,51:60],na.rm = T)
test.all.pa$mean_delta_RUD <- test.all.pa$mean_RUD_ALL -test.all.pa$mean_RUD_CTRL

lengthened_genes <- test.all.pa[mean_delta_RUD>0,]$PA.ID
shortened_genes <- test.all.pa[mean_delta_RUD <= 0,]$PA.ID

L_genes <- unique(do.call(c,lapply(strsplit(lengthened_genes,'\\|'), function(x)x[1])))
S_genes <- unique(do.call(c,lapply(strsplit(shortened_genes,'\\|'), function(x)x[1])))

GO_genes <- unique(do.call(c,lapply(strsplit(all.filtered.id,'\\|'), function(x)x[1])))

genes.names <- unique(do.call(c,lapply(strsplit(all.filtered.id,'\\|'), function(x)x[2])))



library(org.Hs.eg.db)
library(clusterProfiler)
library(Rgraphviz)
PA.GO1 <- enrichGO(gene = GO_genes,
                  OrgDb  = org.Hs.eg.db, 
                  keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
dotplot(PA.GO,showCategory = 20,font.size = 10)
symbols <- genes.names
KEGG_genes<-unique(mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL'))
PA.GO2 <- enrichGO(gene = KEGG_genes,
                  OrgDb  = org.Hs.eg.db, 
                 # keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)


PA.KEGG <- enrichKEGG(KEGG_genes,organism='hsa')

pdf(paste('GO_entrizfiltered.pdf'),height = 8,width = 10)
dotplot(PA.GO2,showCategory = 20,font.size = 10)
plotGOgraph(PA.GO2,firstSigNodes = 15, useInfo = "all")
dev.off()



test.genes <- PA.GO1[PA.GO1$Description == 'lymphocyte differentiation',]$geneID
test.genes <- unlist(strsplit(test.genes,'\\/'))


region <- merge(files[PA.ID%in%all.filtered.id,],ref,by='ID',all.x=T)

region$gene <- do.call(c,lapply(strsplit(region$ID,'\\|'), function(x)x[1]))

region.for.sa <- region[gene%in%test.genes,]


regions <- paste0(region.for.sa$chr,':',region.for.sa$start,':',
                  region.for.sa$end,'.',region.for.sa$gene,'.',region.for.sa$PA)


write.table(regions,file='../sashimiplots/test.region',col.names=F,row.names=F,quote=F)

#'"ENSG00000141968|VAV1|19|+|6857291"'
boxplot(as.numeric(region.for.sa[31,42:51]),as.numeric(region.for.sa[31,52:61]))
#"ENSG00000160685|ZBTB7B|1|+|155016485"
boxplot(as.numeric(region.for.sa[38,42:51]),as.numeric(region.for.sa[38,52:61]))
#"ENSG00000105438|KDELR1|19|-"
#"ENSG00000162739|SLAMF6|1|-"

library(ggplot2)
library(ggsci)
colors = pal_nejm()(7)[3:4]

df <- melt(region.for.sa[38,42:61])
df <- melt(region.for.sa[16,42:61])
df <- melt(region.for.sa[20,42:61])
df <- melt(region.for.sa[39,42:61])
df <- melt(region.for.sa[51,42:61])

df$value[df$value>1] <- NA
df$group <- c(rep('ALL',10),rep('CTRL',10))

p <- ggplot(df,aes(x=group,y=value))+geom_boxplot(fill=colors)+
    geom_signif(comparisons = list(c("ALL", "CTRL")),
              vjust = -0.5,textsize=4,
              test="wilcox.test")+
  labs(y = 'RUD', x = '')+
  labs(title = "SLAMF6")+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold.italic', size = 16),
        legend.position = 'top', legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

                                                                                                    
ggsave(file='SLAMF6.pdf',height =6,width = 6)


df <- melt(region.for.sa[31,42:61])

df$value[df$value>1] <- NA
df$group <- c(rep('ALL',10),rep('CTRL',10))

p <- ggplot(df,aes(x=group,y=value))+geom_boxplot(fill=colors)+
  geom_signif(comparisons = list(c("ALL", "CTRL")),
              vjust = -0.5,textsize=4,
              test="wilcox.test")+
  labs(y = 'RUD', x = '')+
  labs(title = "VAV1")+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold.italic', size = 16),
        legend.position = 'top', legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave(file='VAV1.pdf',height =6,width = 6)



####gene exoression of C/P factors
factors <- c('CPSF1','CPSF2','CPSF3','CPSF4','FIP1L1',
             'CSTF1','CSTF3','CPSF6','CSTF2','CSTF2T','NUDT21',
             'PCF11','SYMPK','PAPOLG','WDR33')
symbols <- factors
genes<-mapIds(org.Hs.eg.db, symbols,'ENSEMBL','SYMBOL')
GeneExprs <- readRDS('../GeneExp/GeneExprs.rds')

use.data <- GeneExprs$counts[genes,]

bam.names <- do.call(c,lapply(strsplit(sample.info$Bam,'\\/'), function(x)x[9]))
bam.names <- gsub('_','.',bam.names)
identical(bam.names,colnames(use.data))



colnames(use.data) <- sample.info$sample
s1 <- log2(use.data[,1]/use.data[,11])
s2 <- log2(use.data[,2]/use.data[,12])
s3 <- log2(use.data[,3]/use.data[,13])
s4 <- log2(use.data[,4]/use.data[,14])
s5 <- log2(use.data[,5]/use.data[,15])
s6 <- log2(use.data[,6]/use.data[,16])
s7 <- log2(use.data[,7]/use.data[,17])
s8 <- log2(use.data[,8]/use.data[,18])
s9 <- log2(use.data[,9]/use.data[,19])
s10 <- log2(use.data[,10]/use.data[,20])




dt <- cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)

rownames(dt) <- factors
pheatmap(dt,cluster_cols = T,cluster_rows = T)

library("DESeq2")
use.data <- GeneExprs$counts
colnames(use.data) <- sample.info$sample
rownames(sample.info) <- sample.info$sample
dds <- DESeqDataSetFromMatrix(countData = use.data,
                              colData = sample.info,
                              design = ~ Group)

dds$Group <- factor(dds$Group, levels = c("Control","ALL"))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)

saveRDS(res,file='GeneExprs.res.rds')
res <- readRDS('GeneExprs.res.rds')
res.use <- res[genes,]
rownames(res.use) <- factors
res.use[res.use$padj >=0.05,]$log2FoldChange <- 0
data <- res.use$log2FoldChange
names(data) <- factors
pdf('CP factors gene expression.pdf',8,6)
barplot(sort(data,decreasing=T),ylab='Log2 fold change(ALL/CTRL)',
        col = pal_nejm()(7),cex.names = 0.5,
        main= paste0('Gene expression fold-change of polyadenylation factors','\n','(FDR<0.05)'))
dev.off()

library(dplyr)
##L_genes and S_genes
intersect_genes <- intersect(L_genes,S_genes)
L_genes <- L_genes[!L_genes%in%intersect_genes]
L_gene_exprs <- res[L_genes,]
L_gene_exprs <- na.omit(L_gene_exprs)
L_gene_exprs <- data.frame(Log2FC=L_gene_exprs$log2FoldChange,
                           padj=as.numeric(L_gene_exprs$padj))

L_gene_exprs<- L_gene_exprs%>%mutate(sig=ifelse(padj <=0.05 ,
                                                ifelse(Log2FC <0,'Down','Up')
                                                ,'Not_DE'))
p1 <- ggplot(L_gene_exprs,aes(x=Log2FC,y=-log10(padj),col=sig))+
  geom_point()+scale_color_manual(values=c("#56B4E9","#999999", "#E69F00"))+
  labs(title = 'Lengthening')

p1
S_genes <- S_genes[!S_genes%in%intersect_genes]
S_gene_exprs <- res[S_genes,]
S_gene_exprs <- na.omit(S_gene_exprs)
S_gene_exprs <- data.frame(Log2FC=S_gene_exprs$log2FoldChange,
                           padj=as.numeric(S_gene_exprs$padj))

S_gene_exprs<- S_gene_exprs%>%mutate(sig=ifelse(padj <=0.05 ,
                                                ifelse(Log2FC <0 ,'Down','Up'),
                                                       'Not_DE'))
p2 <- ggplot(S_gene_exprs,aes(x=Log2FC,y=-log10(padj),col=sig))+
  geom_point()+scale_color_manual(values=c("#56B4E9","#999999", "#E69F00"))+
  labs(title = 'Shortening')


library(cowplot)
plot_grid(p1,p2)
ggsave(file='DE_genes1.pdf',width = 10,height = 4)
df_L1 <- as.data.frame(table(L_gene_exprs$sig))[c(1,3),]
p11 <- ggplot(df_L1,aes(Var1,Freq,fill=Var1))+geom_col()+
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+labs(title = 'Lengthening')

df_S2 <- as.data.frame(table(S_gene_exprs$sig))[c(1,3),]
p22 <- ggplot(df_S2,aes(Var1,Freq,fill=Var1))+geom_col()+
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+labs(title = 'Shortening')
plot_grid(p11,p22)
ggsave(file='DE_genes.pdf',width = 10,height = 4)




#Nudt21 in each sample

NUDT21<-mapIds(org.Hs.eg.db, symbols,'ENSEMBL','SYMBOL')


vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("sample"))





