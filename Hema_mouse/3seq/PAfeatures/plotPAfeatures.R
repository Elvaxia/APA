setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/3seq/PAfeatures')
library(data.table)
library(GenomicFeatures)
library(dplyr)
all.pa <- fread('annotated.PA.txt')
#as pseudogene has polyA tail, we exclude those genes
all.pa <- all.pa[!grep('pseudogene',all.pa$gene_biotype),]
all.pa <- all.pa[!duplicated(all.pa$PA),]#some PA site on different genes
all.pa <- makeGRangesListFromDataFrame(all.pa,keep.extra.columns=T)


##overlapped with annotated PAs
annotated.pa <- fread('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/PAS/annoPA.bed',
                      header = F)
colnames(annotated.pa) <- c('chr','start','end','source','score','strand')
annotated.pa <- annotated.pa%>%mutate(start=start-24,
                                      end=end+24)
annotated.pa <- makeGRangesFromDataFrame(annotated.pa,keep.extra.columns = T)

overlapped <- findOverlaps(all.pa,annotated.pa)

all.pa <- as.data.table(all.pa)

all.pa$Annotated <- 'Novel'
all.pa[unique(overlapped@from),]$Annotated <- 'Annotated'
#as the most two gene_biotype is protein_coding(81.35%,14958),
#lincRNA
#so group other types
genetype <- all.pa$gene_biotype
genetype[is.na(genetype)] <- 'intergenic'
genetype[!genetype%in%c('lincRNA','protein_coding','intergenic')] <- 'other'
all.pa$gene_biotype <- genetype

write.table(all.pa[,-c(1:2)],file='all.pa.annotated.uniqued.txt',
            col.names = T,row.names = F,quote=F,sep='\t')

write.table(all.pa[,-c(1:2)],file='all.pa.annotated.uniqued.tmp',
            col.names = T,row.names = F,quote=F,sep='\t')
#APA_gene protein coding
apa_gene <- all.pa[gene_biotype=='protein_coding',]%>%
  group_by(gene_id)%>%count()%>%filter(n>1)%>%pull(gene_id)

table(all.pa[gene_id%in%apa_gene,]$PA.type)

#plot overlapped pie chart and PA.type,gene,type distribution

library(ggplot2)
library(ggsci)
library(ggrepel)

APA.colors = pal_nejm()(8)
p.pa.type <- melt(all.pa[,c('PA.type','Annotated')])%>%table%>%as.data.frame%>%
            ggplot(aes(PA.type,Freq,fill=Annotated))+
            geom_bar(position = 'dodge', stat='identity')+
            scale_fill_manual(values = APA.colors)+
            geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)+
            labs(y = 'number of PAs')+
            theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 16),
                  legend.position = 'top', legend.title = element_blank(),
                  legend.text = element_text(size=15),
                  axis.text.x = element_text(size=12),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=15),
                  axis.text.y  = element_text(size=15),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
  
pdf('./PAfeatures.Fig/PA.type.pdf',8,6)
p.pa.type
dev.off()



count.data <- as.data.frame(table(all.pa$gene_biotype)) %>%
  arrange(desc(Var1))%>%
  mutate(pro=round(Freq/sum(Freq),2))%>%
  mutate(lab.ypos=cumsum(pro)-0.5*pro)

genetype <- ggplot(count.data, aes(x = "", y = pro, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 1)+
  scale_fill_manual(values = APA.colors[3:6]) +
geom_text_repel(aes(x = 1.4, y = lab.ypos, label = paste(Var1,Freq)),
                size=8,
                                 nudge_x = 0.25,
                                 segment.size = .1,
                                 show.legend = FALSE)+
 theme_void()+ theme(legend.position = "none")

pdf('./PAfeatures.Fig/PA.genetype.pdf',6,6)
genetype
dev.off()



#basedistri
bd <- fread('../PAfeatures/Basedistri/test.txt')
pdf('./PAfeatures.Fig/basedistri.pdf',6.5,6)
plot(bd$A,type='l',lwd = 2,col=APA.colors[1],ylim=c(0,0.65),xlab ='Genome position to PA sites',
     xaxt = 'n',ylab='Frequency')
axis(1, at=c(0,25,50,75,100), labels=c(0,25,'PA site',75,100))
lines(bd$C,col=APA.colors[3],lwd = 2)
lines(bd$G,col=APA.colors[2],lwd = 2)
lines(bd$T,col=APA.colors[4],lwd = 2)
legend("topleft", c('A','G','C','U'), cex = 1,
       fill = APA.colors[1:4],bty = "n")
dev.off()


#UTR length per gene
utr.length <- fread('../../bulk_tianjin/src/utr_length.txt')
utr.gene <- fread('../../bulk_tianjin/src/reduced.3utr.txt')
apa_gene <- fread('all.pa.annotated.uniqued.txt')
utr.apa.gene <- apa_gene[PA.type%in%c('3_UTR','3_UTR_extended'),]
utr.apa.gene <- makeGRangesListFromDataFrame(utr.apa.gene,
                                         keep.extra.columns=T,
                                         split.field ='gene_id' )

utr.gene <- makeGRangesListFromDataFrame(utr.gene,
                                         keep.extra.columns=T,
                                         split.field ='gene_id' )








