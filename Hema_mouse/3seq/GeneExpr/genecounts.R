library(Rsubread)
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/quantifyPA/align/bam')
bamfiles <- system('ls *bam',intern = T)
Exprs <- featureCounts(files=bamfiles,
                    annot.ext = '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/Mus_musculus.GRCm38.sorted.gtf.gz',
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id")

saveRDS(Exprs,file='/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/GeneExpr/GeneExprs.rds')
