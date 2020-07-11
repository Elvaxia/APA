#sample.info
library(readxl)
library(dplyr)
tmp <- as.data.frame(read_excel('./src/179648-2 (1).xlsx', 2))[,1:4]
filepath <- '/mnt/raid63/download/encode/KD-RNA-seq/fastq/align/'

tmp <- tmp[tmp$Cell_line == 'K562',]
tmp <- tmp[tmp$Cell_line == 'HepG2',]

colnames(tmp)[3] <- 'KD_exp'
kd_rep1 <- tmp[,c(1,3)] %>% mutate(
    Bam = paste0(filepath, KD_exp,'_rep1.Aligned.sortedByCoord.out.bam'),
    Group = 'KD',
    rep =1)
kd_rep2 <- tmp[,c(1,3)] %>% mutate(
    Bam = paste0(filepath, KD_exp,'_rep2.Aligned.sortedByCoord.out.bam'),
    Group = 'KD',
    rep =2)


ctrl_rep1 <- tmp[,c(1,4)] %>% mutate(
    Bam = paste0(filepath, Control_exp,'_rep1.Aligned.sortedByCoord.out.bam'),
    Group = 'Ctrl',
    rep=1)
ctrl_rep2 <- tmp[,c(1,4)] %>% mutate(
    Bam = paste0(filepath, Control_exp,'_rep2.Aligned.sortedByCoord.out.bam'),
    Group = 'Ctrl',
    rep=2)

all <- rbind(kd_rep1[,-2],kd_rep2[,-2],ctrl_rep1[,-2],ctrl_rep2[,-2])
all$sample <- paste(all$RBP,all$Group,all$rep,sep='_')

seq_dep <- read.table('./src/sample_seq_depth',header=F)
seq_dep$V1 <- paste0(filepath,seq_dep$V1)
all <- merge(all,seq_dep,by.x='Bam', by.y = 'V1',all.x=T,all.y=F)
colnames(all)[6] <- 'Seq_depth'
all <- all[!is.na(all$Seq_depth),]

write.table(all, file='RBP.sample.info.txt', col.names=T, row.names=F,quote=F,sep='\t')