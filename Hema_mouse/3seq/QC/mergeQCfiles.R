#merge QC jason files
setwd('/mnt/raid61/APA/Project/HSC_mm/QC')
library(rjson)
jsons <- system('ls *json',intern=TRUE)
json.files <- list()
for(i in jsons){
	res <- fromJSON(file = i)
	res <- as.data.frame(unlist(res))
	res <- head(res,22)
	colnames(res)[1] <- gsub('.json','',i)
	json.files[[i]] <- res
}
json.files <- do.call(cbind,json.files)
json.files <- t(json.files)

#combine bowtie output'/mnt/raid61/PRJNA257488_2014science/bam_files/bam_qc_files'
bowtie.res <- list()
for(i in rownames(json.files)){
	res <- read.table(paste('/mnt/raid61/PRJNA257488_2014science/bam_files/bam_qc_files',
		paste(i,'fastq.bam_sorted.stat',sep='.'),sep='/'),fill=T)
	res <- data.frame(bowtie_total.reads=res$V1[1],bowtie_mapped.reads=res$V1[5],bowtie_mapped.perc=res$V5[5])
	bowtie.res[[i]] <- res
}
bowtie.res <- do.call(rbind,bowtie.res)


All.stat <- cbind(json.files,bowtie.res)

sample.info	 <- read.table('/mnt/raid61/APA/Project/HSC_mm/sample.info.txt',sep='\t',header=T)
filnal <- merge(All.stat,sample.info,by.x='row.names',by.y='Run')
write.csv(filnal,file='/mnt/raid61/APA/Project/HSC_mm/QC_stat.csv')
