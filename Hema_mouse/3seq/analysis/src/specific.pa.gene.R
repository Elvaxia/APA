#
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/analysis')
load('rawdata/specific.pa.RData')
geneExp <- readRDS('../GeneExpr/GeneExprs.rds')
geneExp <- geneExp$counts
sample.info <- fread('../sample.info.txt')
colnames(geneExp) <- gsub('.sorted.bam','',colnames(geneExp))
sample.info <- sample.info[match(colnames(geneExp),Run),]
colnames(geneExp) <- sample.info$ColName


##
use.pa.anno <- use.pa.anno[!duplicated(use.pa.anno$gene_id),]
GO.res.genes <- c('Bcl6','Bcr','Cd44','Cdk6','Ezh2','Fech','Foxo3','Gata2','Kit','Klf2','Lyar','Mef2c','Myb','Ncor1',
'Rhd','Stat1','Trim10','Wdr48','Akap8','Cct5','Cct6a','Cdc27','Dkc1','Fen1','H2afy','Mki67','Mnat1','Myc','Ogt',
'Slk','Kat6a','Pip4k2a','Bop1','Ddx21','Ncl','Riok3','Rpf2','Rpl35','Sirt7','Wdr43','Hsp90aa1','Ptprc',
'Rpl22','Tmem131l','Vsir','Zfp609','Celf2','Cnot6l','Gemin5','Nop58','Rplp0','Atic','Atp5a1','Atp5f1','Dut','Gart','Nudt4',
'Ola1','Shmt2','Umps','Mtr','Baz1b','Jmjd1c','Kdm7a','Msl1','Ifnar2','Irak1','Samhd1','Calr','Erap1','Gm7030','H2−D1','March8',
'Tapbp','Acvrl1','Eng','Erdr1','Foxo1','Smc1a')

gene_id <- unique(use.pa.anno[use.pa.anno$gene_name%in%GO.res.genes,]$gene_id)
use.gene.exp <- geneExp[gene_id,]
cells <- unique(sample.info$CellType)
bulk.gene <- list()
for(i in cells){
    cols <- sample.info[CellType==i,]$ColName
    tmp <- use.gene.exp[,colnames(use.gene.exp)%in%cols]
    if(is.null(nrow(tmp))){
        res <- tmp
    }else{
        res <- rowMeans(tmp)
    }
    bulk.gene[[i]] <- res   
}

bulk.gene <- do.call(cbind,bulk.gene)

tmp <- use.pa.anno[use.pa.anno$gene_name%in%GO.res.genes,]
tmp <- tmp[!duplicated(tmp$gene_id),]
tmp <- tmp[match(gene_id,tmp$gene_id),]

rownames(bulk.gene) <- tmp$gene_name
library(pheatmap)
pdf('figures/GO.gene.expression.pdf',width=6,height=10)
 pheatmap(bulk.gene,scale='row',show_rownames=T,fontsize = 8)
 dev.off()