#Generate gene expression matrix
library("Rsubread")
library(data.table)
sampleinfo <- fread('../ALL_sample.info.txt',header=T)
bamfiles <- sampleinfo$Bam
Exprs <- featureCounts(files=bamfiles,
                    annot.ext = '../ref/Homo_sapiens.GRCh38.93.gtf.gz',
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id")
saveRDS(Exprs,file='GeneExprs.rds')

