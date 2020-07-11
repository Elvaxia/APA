library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)


query <- GDCquery(
    project = "TARGET-ALL-P1",   # https://portal.gdc.cancer.gov/
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - FPKM",
    access = "open",
)

GDCdownload(query)
data <- GDCprepare(query)


saveRDS(data,file='ALL_P1.TARGET.Rds')
#counts <- assay(data)
#clin <- colData(data)
#TCGAanalyze_survival(clin, "gender")




expression = t(assay(data))
gene_names <- as.character(sapply(colnames(expression), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
expression <- merge(as.data.frame(colData(data)), expression, by = 0)
for (i in colnames(expression)) {
    if (i %in% gene_names) {
        expression[, i] <- ifelse(expression[, i] > median(expression[, i]),'high','low')
    }
}


genes <- c('ENSG00000111605',#CPSF6
            'ENSG00000136536',#MARCH7
            'ENSG00000175643',#RMI2
            )
surv <- TCGAanalyze_survival(
    expression, 
    clusterCol = 'ENSG00000111605',
    risk.table = F,     # 去掉table
    conf.int = F,       # 去掉误差线
    surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'test.sur.pdf'
)
#使用expression里days_to_last_follow_up进行分组检验，并作图
temp = surv$plot$layers[[6]]$data
max_x = max(surv$plot$data$time) * 0.8
p <- surv$plot + geom_text(aes(x = max_x, y = 0.9), label = i, size = 7)
ggplot2::ggsave(filename = paste0(i, ".pdf"), width = 10, height = 5, plot = p, dpi = 600)
col_names <- c("pvalue",  as.character(temp[, "type"]), "gene")
res <- c(
    as.numeric(gsub("p = ", "", surv$plot$layers[[4]]$aes_params$label)),
    as.numeric(temp[, "x1"]),
    i
)

#由于p2和p3无法prepare成data，只好直接从网站下载手动合并clinical 和 RPKM


prepareData <- function(path){
    RPKM.file.id <- read.table(file.path(path,'analyte.tsv'),sep='\t',fill=T,header=T)
    sample.sheet <- read.table(file.path(path,'gdc_sample_sheet.2020-04-29.tsv'),sep='\t',fill=T,header=T)
    RPKMs <- list()
    for(i in 1:nrow(sample.sheet)){
        files.path <- paste(sample.sheet[i,]$File.ID,sample.sheet[i,]$File.Name,sep='/')
        file <- fread(file.path(path,files.path))
        colnames(file) <- c('gene_id',as.character(sample.sheet$Sample.ID)[i])
        RPKMs[[i]] <- file
    }
    all.RPKM <- Reduce(function(x,y)merge(x,y,by='gene_id'),RPKMs)
    all.RPKM <- t(all.RPKM)
    colnames(all.RPKM) <- all.RPKM[1,]
    all.RPKM <- as.data.frame(all.RPKM[-1,])
    names <- rownames(all.RPKM)
    new.names <- as.character(RPKM.file.id[match(names,RPKM.file.id$sample_submitter_id),]$case_id)
    all.RPKM$case_id <- new.names
    clinical <- read.table(file.path(path,'clinical.tsv'),sep='\t',header=T)
    data <- merge(clinical,all.RPKM,by='case_id')
    return(data)
}

p3_data <- prepareData('TARGET/ALL_p3_RPKM')
p2_data <- prepareData('TARGET/ALL_p2_RPKM')
p1_data <- prepareData('TARGET/ALL_p1_RPKM')
save(p1_data,p2_data,p3_data,file='TARGET_ALL.RData')


load('TARGET_ALL.RData')
colnames(p1_data)[182:60664] <- gsub('\\..+','',colnames(p1_data)[182:60664])

colnames(p2_data)[182:60664] <- gsub('\\..+','',colnames(p2_data)[182:60664])
colnames(p3_data)[182:60664] <- gsub('\\..+','',colnames(p3_data)[182:60664])

#plot
genes <- c('ENSG00000111605', #CPSF6
            'ENSG00000136536', #MARCH7
            'ENSG00000175643', #RMI2
            'ENSG00000203709', #MIR29B2CHG (MIR29B2 and MIR29C host gene)
            'ENSG00000197548', #ATG7
            'ENSG00000207607', #MIR200
            'ENSG00000284214', #MIR29C
            'ENSG00000116001', #TIA1
            'ENSG00000100320' #RBFOX2
            )

p1_data$stage <- 'p1'
p2_data$stage <- 'p2'
p3_data$stage <- 'p3'



use.data <- rbind(p1_data,p2_data,p3_data,stringsAsFactors=F)

use.data$days_to_death <- as.character(use.data$days_to_death)
for (i in colnames(use.data)) {
    if (i %in% genes) {
        use.data[,i]<- as.numeric(use.data[,i])
        use.data[, i] <- ifelse(use.data[, i] > median(use.data[, i]),'high','low')
    }
}

# use.data1 <- p1_data
# for (i in colnames(use.data1)) {
#     if (i %in% genes) {
#         use.data1[,i]<- as.numeric(use.data1[,i])
#         use.data1[, i] <- ifelse(use.data1[, i] > median(use.data1[, i]),'high','low')
#     }
# }

# use.data2 <- p2_data
# for (i in colnames(use.data2)) {
#     if (i %in% genes) {
#         use.data2[,i]<- as.numeric(use.data2[,i])
#         use.data2[, i] <- ifelse(use.data2[, i] > median(use.data2[, i]),'high','low')
#     }
# }

# use.data3 <- p3_data
# for (i in colnames(use.data3)) {
#     if (i %in% genes) {
#         use.data3[,i]<- as.numeric(use.data3[,i])
#         use.data3[, i] <- ifelse(use.data3[, i] > median(use.data3[, i]),'high','low')
#     }
# }

# use.data <- rbind(use.data1,use.data2,use.data3,stringsAsFactors=F)
# use.data$days_to_death <- as.character(use.data$days_to_death)

use.data$use.col <- "separate"

use.data[use.data$ENSG00000203709=='high'&use.data$ENSG00000136536=='low',]$use.col <- "MARCH7+miR-29"


#
use.data.t <- use.data[use.data$primary_diagnosis=='Precursor B-cell lymphoblastic leukemia'|use.data$primary_diagnosis=='B lymphoblastic leukemia/lymphoma, NOS',]

surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'ENSG00000203709',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'test.sur.pdf'
)
surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'ENSG00000136536',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'march7.sur.pdf'
)

surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'use.col',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'both.sur.pdf'
)


surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'ENSG00000175643',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'march7.sur.pdf'
)


surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'ENSG00000116001',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'TIA1.sur.pdf'
)

surv <- TCGAanalyze_survival(
    use.data.t, 
    clusterCol = 'ENSG00000111605',
    risk.table = F,     # 去掉table
    conf.int = F ,     # 去掉误差线
    #surv.median.line = "hv",   # 添加两条median survival days线
    filename = 'CPSF6.sur.pdf'
)