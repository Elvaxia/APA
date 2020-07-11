library(Rsubread)
library(data.table)
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin')
sample.info <- fread('Hema_sampleinfo.txt',header=T)
bamfiles <- sample.info$Bam
Exprs <- featureCounts(files=bamfiles,
                    annot.ext = '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz',
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    isPairedEnd= TRUE,
                    nthreads = 10)

saveRDS(Exprs,file='/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/GeneExprs/GeneExprs.rds')

##
library(DESeq2)
library(data.table)
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin')
Exprs <- readRDS('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/GeneExprs/GeneExprs.rds')
Exprs.counts <- Exprs$counts
colnames(Exprs.counts) <- gsub('.Aligned.sortedByCoord.out.bam','',colnames(Exprs.counts))
colnames(Exprs.counts) <- gsub('[.]','_',colnames(Exprs.counts))
colnames(Exprs.counts)[1:4] <- c('LT-HSC_1','LT-HSC_2','ST-HSC_1','ST-HSC_2')
library(rtracklayer)
gtf <- as.data.table(import('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz'))
gene.name <- gtf[,c('gene_id','gene_name')]
gene.name <- unique(gene.name)
apa.factor <- read.table('./src/apa.factors')
apa.gene <- merge(apa.factor,gene.name,by.x='V1',by.y='gene_name',all.y=F)
colnames(apa.gene)[1] <- 'gene_name'


sample.info <- fread('Hema_sampleinfo.txt',header=T)
rownames(sample.info) <- sample.info$sample
dds <- DESeqDataSetFromMatrix(countData = Exprs.counts,
                              colData = sample.info,
                              design = ~ CellType)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$Celltype <- relevel(dds$CellType, ref = 'LT-HSC')
dds <- DESeq(dds)

res.list <- list()
for(i in unique(sample.info$CellType)[-1]){
    res.list[[i]] <- results(dds,contrast =c('CellType', i,'LT-HSC'))
}
saveRDS(res.list, file='./Analysis/rawdata/DESeq.res.rds')

#extract apa factor genes
genes <- apa.gene$gene_id
apa.data <- as.data.frame(do.call(cbind,lapply(res.list, function(x){
    te <- as.data.frmae(x[genes,c(2,6)]
    te <- te[te$padj <0.1,]
})))

filtered.res <- list()
for(i in 1:16){
    te <- as.data.frame(res.list[[i]][genes,c(2,6)])
    te <- te[te$padj <0.1,]
    te$gene <- rownames(te)
    te <- te[,c(1,3)]
    filtered.res[[i]] <- te
}
filtered.merged <- Reduce(function(x,y){
    merge(x,y,by='gene',all=T)
},filtered.res)

colnames(filtered.merged) <- c('gene',names(res.list))
gene.names <- apa.gene[match(filtered.merged$gene, apa.gene$gene_id),]$gene_name
apa.data <- filtered.merged[,c(-1,-2,-3,-17)]
rownames(apa.data) <- gene.names
nas <- apply(apa.data, 1, function(x)sum(is.na(x)))
apa.data <- apa.data[nas<10,]

apa.data[is.na(apa.data)]<- 0
library(pheatmap)
pdf('./Analysis/figures/apa.factor.foldchange.pdf',width=10,height=6)
pheatmap(t(apa.data),border_color = NA,fontsize_row = 12,fontsize_col =5.5,cellwidth =6,cellheight=15,cluster_col=T,color = rev(c('#8e0152','#c51b7d','#de77ae','#f1b6da','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221')))
dev.off()

