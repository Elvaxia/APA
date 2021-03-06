#use hypergemometric test find cell type specific genes
library(here)
library(data.table)
library("factoextra")
library("FactoMineR")
library(dplyr)
library(dendextend)
library(ggplot2)

RUCP <- fread(here('/Analysis/rawdata//filtered.RUD.txt'))
rucp <- as.data.frame(RUCP[,130:193])
rownames(rucp) <- paste0(RUCP$gene_id.x,":",RUCP$PA.ID)

sample.info <- fread('hs_hema_sampleinfo.txt')






anno_col <- gsub('_.*','',colnames(rucp))
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
  i.hyper.pa<-i.hyper[i.hyper$adj_pval < 0.01 & i.hyper$Cell_p > .1 & i.hyper$Cell_p > 2*i.hyper$Other_p, ]$id
  df<-data.frame(cells=rep(a,length(i.hyper.pa)),pa=i.hyper.pa)
  return(df)
}
#lengthened
cell.specific <- list()
for(i in unique(anno_col)){
  cell.specific[[i]] <- hyper_func(i,rucp)
}
cell.specific <- cell.specific[1:16]

all.res <- do.call(c,lapply(cell.specific,function(x)as.character(x$pa)))
length(all.res)#2895
rows <- unique(all.res)
length(rows)#1579
plot.data <- rucp[rownames(rucp)%in%rows,]
mean_rud <- list()

for(i in unique(anno_col)){
    col.n <- which(anno_col == i)
  mean_rud[[i]] <- rowMeans(plot.data[,col.n])
}
mean_rud <- do.call(cbind,mean_rud)

colnames(mean_rud) <- unique(anno_col)

mean_rud <- mean_rud[,1:16]
#plot

ge<-t(scale(mean_rud))
res <- hcut(ge,stand = T)
# Visualize
dend <- as.dendrogram(res)

pdf(here('Analysis/figures/celltypespecific.distance.pdf'),height=6,width=3)
dend %>% set("leaves_pch", 15) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", c('#8dd3c7',rep('#ffffb3',3),rep('#fccde5',2),rep('#bebada',2),
  rep('#fb8072',2),
  rep('#80b1d3',2),
  rep('#fdb462',2),
  rep('#33a02c',2))) %>% # node point color
  plot(main = "Distance between cell types",horiz=T)
dev.off()

library(pheatmap)
library(RColorBrewer)
pdf(here('Analysis/figures/celltypespecific.heatmap.pdf'),height=4,width=6)
pheatmap(t(ge), show_rownames  = FALSE, clustering_method='ward.D2',scale='row',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(10))
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


names <- c('HSC','GMP','CMP','MEP','CLP','MPP','DC','Mp','EB','MK','Mo','Neu','CD4','CD8','B','NK')
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
labs(x='',y='Number of specifically detcted APA events',fill='Specifically detected in groups') +
scale_y_continuous(expand = c(0, 0)) +
coord_flip()
dev.off()

save(events.df,file=here('Analysis/rawdata/hs.group.detected.RData'))
