setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP')
library(here)
library(data.table)
library(dplyr)
library(ggplot2)
require("ggrepel")
load(here('Analysis/rawdata/rbp.testres.2.21.RData'))

#save all apa genes of each RBP
extract_apa_genes <- function(rbp.DE.dt){
    rbp.DE.dt[rbp.DE.dt$kd.rud.mean >1,]$kd.rud.mean <- 1
    rbp.DE.dt[rbp.DE.dt$ctrl.rud.mean >1,]$ctrl.rud.mean <- 1
    rbp.DE.dt <- na.omit(rbp.DE.dt)
    rbp.DE.dt$gene <- do.call(c,lapply(strsplit(rbp.DE.dt$id,'\\:'),function(x)x[1]))
    rbp.DE.dt <- rbp.DE.dt %>% mutate(res = ifelse(p_value <0.05,
                                                    ifelse(log2FC >= 0.5,'Lengthened',
                                                    ifelse(log2FC <= -0.5,'Shortened','Not significant')),'Not significant'))
    L.apa.genes <- unique(rbp.DE.dt[rbp.DE.dt$res=='Lengthened',]$gene)
    S.apa.genes <- unique(rbp.DE.dt[rbp.DE.dt$res=='Shortened',]$gene)
    L.apa.genes <- setdiff(L.apa.genes,S.apa.genes)
    S.apa.genes <- setdiff(S.apa.genes,L.apa.genes)
    res.dt <- data.frame(gene = c(L.apa.genes,S.apa.genes),
                        lable = c(rep('Lengthened',length(L.apa.genes)),c(rep('Shortened',length(S.apa.genes)))))
    return(res.dt)
}

hepg2.apa.genes <- lapply(hepg2.de,extract_apa_genes)
k562.apa.genes <- lapply(new.k562.de,extract_apa_genes)

save(hepg2.apa.genes,k562.apa.genes,file=here('Analysis/rawdata/all.apa.genes.2.21.RData'))
#all lengthened genes GO
hepg2.l.genes <- unique(do.call(c,(lapply(hepg2.apa.genes,function(x){
    tmp <- as.character(x[x$lable=='Lengthened',]$gene)
    return(tmp)
})))) 

k562.l.genes <- unique(do.call(c,(lapply(k562.apa.genes,function(x){
    tmp <- as.character(x[x$lable=='Lengthened',]$gene)
    return(tmp)
}))))
overlapped.l.genes <- intersect(hepg2.l.genes,k562.l.genes)

#all shortened genes GO
hepg2.s.genes <- unique(do.call(c,(lapply(hepg2.apa.genes,function(x){
    tmp <- as.character(x[x$lable=='Shortened',]$gene)
    return(tmp)
})))) 

k562.s.genes <- unique(do.call(c,(lapply(k562.apa.genes,function(x){
    tmp <- as.character(x[x$lable=='Shortened',]$gene)
    return(tmp)
}))))
overlapped.s.genes <- intersect(hepg2.s.genes,k562.s.genes)


overlapped.L.genes <- setdiff(overlapped.l.genes,overlapped.s.genes)
overlapped.S.genes <- setdiff(overlapped.s.genes,overlapped.l.genes)
#GO
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(Rgraphviz)
GO_genes <- overlapped.L.genes
#
GO_genes <- c("CD160" , "TYROBP" ,"TRDC"   ,"SIVA1" , "DUSP2" , "RORA"  , "CMC1"  , "YWHAQ",#cluster CD8 sense
                "CD70","CCL4","COTL1","PTGDS","TRDC","IFNG","LGALS1",#cluster CD8 act
                "S100B","IGLC3"#cluster CD8 naive
                )
ensembl.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, keytype = "SYMBOL", column="ENSEMBL")
ensembl.genes <- na.omit(ensembl.genes)
PA.GO.L <- enrichGO(gene = ensembl.genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
pdf('AL.DE.gene.GO.pdf',8,6)
barplot(PA.GO.L)
dev.off()

#

PA.GO.L <- enrichGO(gene = GO_genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
GO_genes <- overlapped.S.genes
PA.GO.S <- enrichGO(gene = GO_genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

symbols <- c(overlapped.L.genes,overlapped.S.genes,dynamic.genes)
entrize_genes<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ENSEMBL')
identical(names(entrize_genes),symbols)
mydf <- data.frame(Entrez = as.numeric(entrize_genes),group = c(rep('Lengthened',length(overlapped.L.genes)),
                                                                rep('Shortened',length(overlapped.S.genes)),
                                                                rep('Dynamic',length(dynamic.genes))))
xx.formula <- compareCluster(Entrez~group, data=mydf,
                              OrgDb='org.Hs.eg.db',
                                  fun='enrichGO' )
dotplot(xx.formula)
dev.off()

dynamic.genes <- intersect(overlapped.l.genes,overlapped.s.genes)
GO_genes <- dynamic.genes
PA.GO.dy <- enrichGO(gene = GO_genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

#CPSF6 target genes GO

cpsf6.genes.hepg2 <- hepg2.apa.genes[['CPSF6']]
cpsf6.genes.k562 <- k562.apa.genes[['CPSF6']]

GO_genes <- as.character(unique(cpsf6.genes.k562$gene))
PA.GO.cpsf6.k562 <- enrichGO(gene = GO_genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
entrize_genes<-mapIds(org.Hs.eg.db, GO_genes, 'ENTREZID', 'ENSEMBL')

PA.KEGG <- enrichKEGG(entrize_genes,organism='hsa')
edox <- setReadable(PA.KEGG, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox)
pdf(here('Analysis/figures/cpsf6.k562.kegg.pdf'),8,8)
p1
dev.off()