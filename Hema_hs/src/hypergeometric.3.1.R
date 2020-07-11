#use hypergemometric test find cell type specific genes
library(here)
library(data.table)
RUCP <- fread(here('/Analysis/rawdata//filtered.RUD.txt'))
rucp <- as.data.frame(RUCP[,130:193])
rownames(rucp) <- paste0(RUCP$gene_id.x,":",RUCP$PA.ID)

sample.info <- fread('hs_hema_sampleinfo.txt')
#remove UVEC mp.ac.mp.inf,and EB,MK ployA RNA
sample.info <- sample.info[!CellType%in%c("UVEC_P","UVEC_R","Mp_act","Mp_inf"),]
#remove EB,MK ployA RNA
total <- sample.info[sample.info$CellType%in%c("EB","MK")&lib!='polyA RNA',]$sample
sample.info <- sample.info[!sample%in%total,]


colnames(rucp) <- gsub('.RUD','',colnames(rucp))
rucp <- rucp[,colnames(rucp)%in%as.character(sample.info$sample)]

anno_col <- gsub('_[1-9].*','',colnames(rucp))
rucp[rucp >1] <- NA
rucp[is.na(rucp)] <- 0



hyper_func<-function(a,DT){
  i.hyper<-t(apply(DT,1,function(x){#change the DT if need
    k=sum(x[which(anno_col==a)]>=0.5 & x[which(anno_col==a)] <=1)
    D=sum(x>=0.15 & x<=1)
    n=length(x)-D
    N=length(which(anno_col==a))
    pval=phyper(k,D,n,N,lower.tail=FALSE)
    if(k==0){
      adj_pval<-pval
    }else{
      adj_pval<-pval*k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval")
    return(enrichment)
  }))

  i.hyper <- as.data.frame(i.hyper)
  i.hyper$Cell_p <- i.hyper$k/i.hyper$N
  i.hyper$Bg_p <- i.hyper$D/(i.hyper$D + i.hyper$n)
  i.hyper$Other_p <- (i.hyper$D - i.hyper$k)/(i.hyper$D + i.hyper$n - i.hyper$N)
  i.hyper$id<-rownames(DT)#change the DT if need
  i.hyper.pa<-i.hyper[i.hyper$adj_pval < 0.01 & i.hyper$Cell_p > .2 & i.hyper$Cell_p > 2*i.hyper$Other_p, ]$id
  df<-data.frame(cells=rep(a,length(i.hyper.pa)),pa=i.hyper.pa)
  return(df)
}
#lengthened
cell.specific <- list()
for(i in unique(anno_col)){
  cell.specific[[i]] <- hyper_func(i,rucp)
}

all.res <- do.call(c,lapply(cell.specific,function(x)as.character(x$pa)))
length(all.res)#3321
rows <- unique(all.res)
length(rows)#1709
plot.data <- rucp[rownames(rucp)%in%rows,]
mean_rud <- list()
for(i in unique(anno_col)){
    col.n <- which(anno_col == i)
  mean_rud[[i]] <- rowMeans(plot.data[,col.n])
}
mean_rud <- do.call(cbind,mean_rud)
colnames(mean_rud) <- unique(anno_col)
#plot
library("factoextra")
library("FactoMineR")
ge<-t(scale(mean_rud))
res <- hcut(ge,stand = T)
# Visualize
dend <- as.dendrogram(res)
library(dplyr)
library(dendextend)
pdf(here('Analysis/figures/celltypespecific.distance.pdf'),height=6,width=3)
dend %>% set("leaves_pch", 15) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", c('#8dd3c7',rep('#ffffb3',2),rep('#fccde5',3),rep('#bebada',2),
  rep('#fb8072',2),
  rep('#80b1d3',2),
  rep('#fdb462',2),
  rep('#33a02c',2))) %>% # node point color
  plot(main = "Distance between cell types",horiz=T)
dev.off()

#detected in myeloid,detected in lymphoid,detetcted in myloied,
#detected in progenitors.detetcted in single celltypes.

all <- do.call(rbind,cell.specific)

dup <- all.res[duplicated(all.res)]
all <- all%>%mutate(enriched=ifelse(pa%in%dup,'multiple','single'))
test.all <- all[all$enriched=='multiple',]
single.all <- all[all$enriched=='single',]
test.events <- test.all%>%group_by(pa)%>%
                mutate(detected = ifelse(cells%in%c('HSC','MPP','CLP','CMP','GMP','MEP'),
                                                        'Progenitors',
                                                        ifelse(cells%in%c('B','CD4','CD8','NK'),'Lymphoid',
                                                            ifelse(cells%in%c('EB','MK'),'Erythroid','Myeloid'))))

single.all$detected <- 'Single'

events.df <- rbind(single.all,as.data.frame(test.events))

#save(events.df,file=here('Analysis/rawdata/specific.apa.df.RData'))

#visualize events
library(ggplot2)

names <- c('HSC','CLP','MPP','CMP','GMP','MEP','DC','Mp','EB','MK','Mo','Neu','CD4','CD8','B','NK')
te <- events.df%>%group_by(cells,detected)%>%count()
pdf(here('Analysis/figures/detcted.groups.pdf'),8,6)
ggplot(te,aes(x=cells,y=n,fill=detected)) + geom_bar(stat="identity",width=0.8) +
scale_fill_manual(labels = c("Erythroid", "Lymhoid",'Myeloid','Progenitors','Single'), 
                    values = c('#beaed4','#7fc97f','#fdc086','#ffff99','#386cb0'))+
scale_x_discrete(limits=names)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        #legend.title=element_blank(),
                        legend.position='right',
                        legend.text=element_text(size=12),
                        axis.title = element_text(size=22),
                        axis.text = element_text(size=12))+
labs(x='',y='Number of detcted APA events',fill='Detected in groups') +
scale_y_continuous(expand = c(0, 0)) +
coord_flip()
dev.off()

save(events.df,file=here('Analysis/rawdata/hs.group.detected.RData'))
