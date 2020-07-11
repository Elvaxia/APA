#
library(here)
library(data.table)
library(rtracklayer)
library(DESeq2)
geneExp <- readRDS(here('GeneExp/GeneExprs.rds'))
geneExp <- geneExp$counts
colnames(geneExp) <- gsub('.Aligned.sortedByCoord.out.bam','',colnames(geneExp))
colnames(geneExp) <- gsub('[.]','_',colnames(geneExp))
all.sample <- fread(here('ALL_sample.info.txt'))

#
RBPs <- fread('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP/RBP.sample.info.txt')
RBPs <- unique(RBPs$RBP)
#
gtf <- as.data.table(import('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/Homo_sapiens.GRCh38.93.sorted.gtf.gz'))
gene.name <- gtf[,c('gene_id','gene_name')]
gene.name <- unique(gene.name)
#
RBPs.gene <- gene.name[match(RBPs,gene.name$gene_name),]

#
rownames(all.sample) <- all.sample$ID
dds <- DESeqDataSetFromMatrix(countData = geneExp,
                              colData = all.sample,
                              design = ~ Group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)


#extract RBPs genes

res <-  results(dds,contrast =c('Group','ALL','Control'))
#ENSG00000167005    NUDT21
#ENSG00000071894     CPSF1
#ENSG00000136536   MARCH7
#ENSG00000175643  RMI2
res.05<- subset(res, padj<0.05)
res.05 <- subset(res.05,log2FoldChange >0.5|log2FoldChange< - 0.5)

select <- c(intersect(RBPs.gene$gene_id,rownames(res.05)),'ENSG00000167005','ENSG00000071894')
RBPs.res <- res.05[select,]
RBPs.res$gene_name <- gene.name[match(select,gene.name$gene_id),]$gene_name
write.table(RBPs.res,file=here('Analysis/rawdata/RBPs.expression.txt'),col.names=T,row.names=T,sep='\t',quote=F)





library(ggplot2)
df <- as.data.frame(RBPs.res)
df <- df%>%mutate(res = ifelse(log2FoldChange>0,'Up','Down'))
name.ordered <- df[order(df$log2FoldChange,decreasing=F),]$gene_name
pdf(here('Analysis/figures/RBPs.expression.pdf'),width=8,height=4)
ggplot(df,aes(x=gene_name,y=log2FoldChange,fill=res))+ geom_col()+
scale_fill_manual(labels = c("Up", "Down"), values = c("#a1d76a","#e9a3c9"))+
scale_x_discrete(limits=rev(name.ordered))+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=0.5, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="top",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text.x =element_text(angle=90),
                        axis.text = element_text(size=8))+
labs(x='',y='Log2(ALL/normal)') +
scale_y_continuous(position = "left",expand = c(0, 0),limits=c(-1,3.1)) 
#coord_flip()
dev.off()

plotCounts(ntd, gene='ENSG00000136536', intgroup="Group")
dev.off()
library("pheatmap")
ntd <- normTransform(dds)
plotMA(res)
dev.off()

df <- counts(dds)[c('ENSG00000136536','ENSG00000175643'),]
df <- assay(ntd)[c('ENSG00000136536','ENSG00000175643'),]
colnames(df) <- all.sample$Group
saveRDS(df, file= './Analysis/rawdata/march7andRMI2geneexpression.rds')
df <- readRDS('./Analysis/rawdata/march7andRMI2geneexpression.rds')

plot.data <- melt(as.matrix(df[1,]))
ggplot(plot.data, aes(x =Var1, y =value, fill = Var1)) + geom_boxplot() +
  stat_compare_means()+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
                  legend.position = 'top', legend.title = element_blank(),
                  legend.text=element_text(size=16),
                  axis.title.x = element_text(size=25),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=25),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
labs( x = '', y='Normalized gene expression')



pheatmap(RBPs.res)
df <- as.data.frame(colData(dds)[,c('ID','Group')])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()
sampleDists <- dist(t(assay(ntd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- ntd$sample
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
p <- plotPCA(ntd,intgroup='sample',returnData=T)
dev.off()
ggplot(p,aes(PC1,PC2,label=sample))+geom_point()+geom_text()
dev.off()
####
library(edgeR)
y <- DGEList(counts=geneExp)
group <- all.sample[match(colnames(geneExp),all.sample$ID),]$Group
y <- DGEList(counts=geneExp, group=group)
keep <- filterByExpr(y)
 y <- y[keep, , keep.lib.sizes=FALSE]
 y <- calcNormFactors(y)
 y <- estimateDisp(y)
 et <- exactTest(y)
 topTags(et)
 res.edgeR <- et$table
 library(dplyr)
 res.RBPs.edegR <- res.edgeR%>%filter(FDR <0.05,abs(logFC) >1.5)