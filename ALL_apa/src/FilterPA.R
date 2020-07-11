setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/ALL_apa')
library(data.table)
library(dplyr)
filenames <- list.files('./ChangePoint/')
files <- list()
for ( i in filenames ) { files[[i]] <- fread(paste0('./ChangePoint/',i) , drop=1)}

files <- do.call(rbind,files)
#filtering out some events
#1.no APA
files <- files[!is.na(PA),]
#2.S_coverage lower than 10,convert RUD to be NA
sample.info <- fread('./ALL_sample.info.txt',header=T)
samples <- sample.info$sample
newfiles <- files

newfiles$PA.ID <- paste(newfiles$gene_id,newfiles$PA,sep = '|')


ref.anno <- fread('./src/reduced.3utr.txt',header = T)
#remove PAS within 150bp distance to the start and end of region
ref.anno.m <- ref.anno%>%group_by(gene_id)%>%mutate(st=min(start),ed=max(end))%>%
  dplyr::select(seqnames,strand,gene_id,st,ed)%>%unique()

PAS <- merge(newfiles[,c('PA','gene_id','PA.ID')],ref.anno.m,by='gene_id',all.x = T)

dist_s <- abs(PAS$PA- PAS$st)
dist_e <- abs(PAS$PA-PAS$ed)

PAS <- PAS[dist_s >=200 &dist_e >=200,]

test.files <- merge(newfiles,PAS,by='PA.ID',all.y = T)
test.files$PA.ID <- paste0(test.files$seqnames,':',test.files$PA.x,':',test.files$strand)
#3.mean RUD of all samples greater than  1.
#4.coverage of cUTR should be over 20
ruddt <- test.files[,42:61]

mean_RUD <- rowMeans(ruddt,na.rm = T)
sd_RUD <- apply(ruddt,1,function(x)sd(x,na.rm=T))
median_RUD <- apply(ruddt,1,function(x)median(x,na.rm=T))
N.NAs <- apply(ruddt,1,function(x)sum(is.na(x)))
S.cov <- apply(test.files[,22:41], 1, median)
S.cov_na <- apply(test.files[,22:41],1,function(x){
  sum(x<=10)
})
merged.genes <- apply(ruddt,1,function(x)sum(x>1.1,na.rm=T))

#6444


filtered <- test.files[S.cov>=20&N.NAs <=18&
                         merged.genes <=2,]
                
filtered$id <- paste(filtered$gene_id.x,
                     filtered$seqnames,
                     filtered$strand,
                     filtered$st,
                     filtered$PA.y,
                     filtered$ed,
                     sep = ':'
)
write.table(filtered,file=here('/Analysis/rawdata/filtered.RUD.txt'),col.names=T,row.names=F,quote=F)

#========================================================================
sashimi.plot(id=filtered$id[69])


#DE apa statistics=======================================================
filtered$mean_RUD_ALL <- rowMeans(filtered[,42:51],na.rm = T)
filtered$mean_RUD_CTRL <- rowMeans(filtered[,52:61],na.rm = T)
filtered$mean_delta_RUD <- filtered$mean_RUD_ALL -filtered$mean_RUD_CTRL
filtered$log2FC <- log2(filtered$mean_RUD_ALL/filtered$mean_RUD_CTRL)
t.test.mtx <- filtered[,42:61]
    t.test.mtx[t.test.mtx > 1] <- NA
    suppressWarnings(p_value <- apply(t.test.mtx,1,function(x){
      if(NA%in%x){
        p.v <- NA
      }else{
        p.v <- t.test(x[1:10],x[11:20])$p.value
      }
      return(p.v)
    }))
filtered$p.value <- p_value
lengthened_genes <- filtered[log2FC>= 0.5 & p.value<0.05,]$id #81
shortened_genes <- filtered[log2FC <= -0.5 &p.value<0.05,]$id #334

L_genes <- unique(do.call(c,lapply(strsplit(lengthened_genes,'\\:'), function(x)x[1])))#75

S_genes <- unique(do.call(c,lapply(strsplit(shortened_genes,'\\:'), function(x)x[1]))) #303
pdf(here('Analysis/figures/L.s.statistic.pdf'),height=6,width=8)
dt.use <- filtered
dt.use <- dt.use %>% mutate(res = ifelse(p.value <0.05,
                                                ifelse(log2FC >= 0.5,'Lengthened',
                                                ifelse(log2FC <= -0.5,'Shortened','Not significant')),'Not significant'))
dt.use <- na.omit(dt.use)
ggplot(dt.use,aes(x= log2FC,y= -log10(p.value),fill=factor(res))) + geom_point(aes(col=factor(res))) +
theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="top",
                        legend.text=element_text(size=16),
                        axis.title = element_text(size=25),
                        axis.text = element_text(size=12))

ggplot(dt.use[dt.use$res!='Not significant',],aes(x=mean_RUD_CTRL,y=mean_RUD_ALL))+
  geom_point(size=1)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=22),
                        axis.text = element_text(size=12))+
                        xlim(c(0,1))+ylim(c(0,1))+labs(x= 'RUCP of control', y='RUCP of ALL')+geom_abline(slope=1, intercept=0) +
 annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(dt.use$res)[1]),color='red',size=6) +
  annotate("text", x = 0.8, y=0.1, label = paste0('Shortened',' ',table(dt.use$res)[3]),color ='blue',size=6)
dev.off()

pdf('hexbin.pdf',8,6)
ggplot(dt.use[dt.use$res!='Not significant',],aes(x=mean_RUD_CTRL,y=mean_RUD_ALL))+
  # define the color of the outline of the hexagons with color=c()
  # using c(#"809FFF") allows for the usage of hexadecimal color codes
  stat_binhex(color=c("#D7DADB")) +
  # set the graph theme to classic, provides a white background and no grid lines
  # Change font size to 18 by using base_size = 18
  theme_classic(base_size=20) +
  # Apply lables to the graph for x and y
  # change the gradient fill to range from grey to Red
  scale_fill_gradient2(low="grey80",high="red")+
  #scale_x_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1))+
  #scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8,1)) +
  geom_abline(slope=1, intercept=0,linetype=3) +
  labs(x= 'RUCP of control', y='RUCP of ALL')+geom_abline(slope=1, intercept=0) +
  annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(dt.use$res)[1]),color='red',size=5.5) +
  annotate("text", x = 0.6, y=0.1, label = paste0('Shortened',' ',table(dt.use$res)[3]),color ='blue',size=5.5)
dev.off()
#################=================
#GO
#GO
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(Rgraphviz)



PA.GO.S<- enrichGO(gene = c(L_genes,S_genes),
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
pdf(here('Analysis/figures/GO.pdf'),height=3,width=8)
barplot(PA.GO.S)  
dev.off()  
PA.GO.S[3,'geneID']
# "ENSG00000071655/ENSG00000092871/ENSG00000101266/ENSG00000102054/ENSG00000132383/ENSG00000133119/ENSG00000135679/ENSG00000136536/ENSG00000163564/ENSG00000163781/ENSG00000175643/
#ENSG00000177374/ENSG00000181090"

#RMI2 ENSG00000175643:16:+:11350654:11351387:11351762
#MARCH7 ENSG00000136536:2:+:159767343:159768515:159771027
#HIC1 ENSG00000177374:17:+:2058836:2061141:2063241


pdf(here('Analysis/figures/MARCH7.pdf'),2.5,4)
df <- filtered[id=='ENSG00000136536:2:+:159767343:159768515:159771027',42:61]
boxplot(as.numeric(df[,1:10]),as.numeric(df[,11:20]),col=c("#5caf90","#f7dc56"),xaxt='n',ylab='RUCP',cex.lab=1,cex.aixs=1)
text(1.5, 0.68, format(filtered[id=='ENSG00000136536:2:+:159767343:159768515:159771027',]$p.value, scientific = TRUE,digits = 2),
cex = 0.8)
axis(1,1:2,c('ALL','Normal'),cex.axis=1)
dev.off()

pdf(here('Analysis/figures/RMI2.pdf'),2.5,4)
df <- filtered[id=='ENSG00000175643:16:+:11350654:11351387:11351762',42:61]
boxplot(as.numeric(df[,1:10]),as.numeric(df[,11:20]),col=c("#5caf90","#f7dc56"),xaxt='n',ylab='RUCP',cex.lab=1,cex.aixs=1)
text(1.5, 0.9, format(filtered[id=='ENSG00000175643:16:+:11350654:11351387:11351762',]$p.value, scientific = TRUE,digits = 2),
cex = 0.8)
axis(1,1:2,c('ALL','Normal'),cex.axis=1)
dev.off()




PA.GO.L<- enrichGO(gene = L_genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
barplot(PA.GO.L)   
dev.off()

symbols <- c(S_genes,L_genes)
entrize_genes<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ENSEMBL')
PA.KEGG <- enrichKEGG(entrize_genes,organism='hsa')
edox <- setReadable(PA.KEGG, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox)
dev.off()
dev.off()              
#============================
#gene expression and s——_
library(rtracklayer)
library(DESeq2)
library(here)
library(data.table)
filtered <- fread(here('/Analysis/rawdata/filtered.RUD.txt'))
geneExp <- readRDS(here('GeneExp/GeneExprs.rds'))
geneExp <- geneExp$counts
colnames(geneExp) <- gsub('.Aligned.sortedByCoord.out.bam','',colnames(geneExp))
colnames(geneExp) <- gsub('[.]','_',colnames(geneExp))
all.sample <- fread(here('ALL_sample.info.txt'))
rownames(all.sample) <- all.sample$ID
dds <- DESeqDataSetFromMatrix(countData = geneExp,
                              colData = all.sample,
                              design = ~ Group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <-  results(dds,contrast =c('Group','ALL','Control'))


ntd <-  normTransform(dds)
saveRDS(ntd,file=here('GeneExp/ntd.rds'))
assay(ntd)['ENSG00000136536',]->march7
assay(ntd)['ENSG00000203709',]->mir29
cor.test(march7,mir29)
#p-value = 0.0003565
#cor -0.7187378 
plot(march7,mir29)



use.genes <- intersect(S_genes,rownames(assay(ntd)))
use.gene<- 'ENSG00000175643'
S.genes.exp <- as.data.frame(assay(ntd))[use.genes,]
identical(colnames(S.genes.exp),sample.info$ID)
colnames(S.genes.exp) <- sample.info$sample

S.pdui <- filtered[log2FC<= -0.5 & p.value<0.05,42:63]

S.pdui <- S.pdui[!duplicated(S.pdui$gene_id.x),]
S.pdui <- as.data.frame(S.pdui[match(use.genes,S.pdui$gene_id.x),])

identical(as.character(rownames(S.genes.exp)),as.character(S.pdui$gene_id.x))
S.genes.exp <- S.genes.exp[S.pdui$gene_id.x,]
#L

use.genes <- intersect(L_genes,rownames(assay(ntd)))
L.genes.exp <- as.data.frame(assay(ntd))[use.genes,]
identical(colnames(L.genes.exp),sample.info$ID)
colnames(L.genes.exp) <- sample.info$sample

L.pdui <- filtered[log2FC>= 0.5 & p.value<0.05,42:63]


L.pdui <- L.pdui[!duplicated(L.pdui$gene_id.x),]
L.pdui <- as.data.frame(L.pdui[match(use.genes,L.pdui$gene_id.x),])

identical(as.character(rownames(L.genes.exp)),as.character(L.pdui$gene_id.x))
L.genes.exp <- L.genes.exp[L.pdui$gene_id.x,]

#

pdf(here('Analysis/figures/geneexp.rucp.pdf'),4,4)
df <- data.frame(PDUI=c(S.pdui[,1],L.pdui[,1]),Gene=c(S.genes.exp[,1],L.genes.exp[,1]))
df <- df[df$Gene!=0,]
corres <- cor.test(df[,1],df[,2])
plot(df,xlab='RUCP',ylab='Normalized gene expression',col='red',pch=19)
text(0.6, 4, paste0('r=',format(corres$estimate, digits = 2)),cex = 0.8)
text(0.6, 3, paste0('p.value=',format(corres$p.value, scientific = TRUE,digits = 2)),cex = 0.8)
dev.off()


library(ggplot2)
library(ggpubr)
library(LSD)
library(MASS)



df <- data.frame(PDUI=c(S.pdui[,1],L.pdui[,1]),Gene=c(S.genes.exp[,1],L.genes.exp[,1]))
df <- df[df$Gene!=0,]
df <- mk_PointDens(df)

  ggplot(data = df, aes(x = PDUI, y = Gene, color = pointdens)) +
    geom_point() +
    coord_cartesian(xlim = c(0, 1), ylim = c(4, 16))+
    scale_colour_gradientn(colours = rev(rainbow(5)))+
    geom_smooth(method = "lm", lwd = 0.5) +
    stat_cor(data = df, method = "pearson",
             size = 6, col = "black",
             label.y = 5, label.x=0.4)+
    labs(subtitle = paste("N =", nrow(df)),
         x = 'RUCP',
         y = 'Normalized gene expression')+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 22),
          axis.line = element_line(),
          axis.text = element_text(size = 16), legend.position = "right")

  ggsave('/mnt/raid61/Personal_data/xiaoxia/APA/Project/ALL_apa/cor.gene.rucp.pdf')
mk_PointDens <-
  function(df,
           x,
           y) {
    x <- df[, x]
    y <- df[, y]
    
    keep <- (!is.na(x)) & (!is.na(y))
    
    x <- x[keep]
    y <- y[keep]
    df <- df[keep, ]
    
    dens <- kde2d(x, y)
    
    # create a new data frame of that 2d density grid
    # (needs checking that I haven't stuffed up the order here of z?)
    gr <- data.frame(with(dens, expand.grid(x, y)), as.vector(dens$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr ~ xgr*ygr, data = gr)
    
    # Apply the model to the original data to estimate density at that point
    df$pointdens <- predict(mod, newdata = data.frame(xgr = x, ygr = y))
    
    return(df)
  }

QuadStac(df, 'PDUI','Gene')
