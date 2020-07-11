library(here)
library(circlize)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
load('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP/Analysis/rawdata/all.apa.genes.2.21.RData')
load(here('Analysis/rawdata/hs.group.detected.RData'))


all.hs.types <- list()
for(i in c(as.character(unique(events.df$cells)),unique(events.df$detected))){
    tmp <- unique(events.df[events.df$cells==i|events.df$detected==i,]$pa)
    tmp <- unique(gsub(':.*','',tmp))
    all.hs.types[[i]] <- tmp
}
all.hs.types <- all.hs.types[-which(names(all.hs.types)=='Single')]

#cluster RBPs by overlapped percentage of all.hs.apa
enrich.lineage <- function(L){
    rbp.n <- do.call(c,lapply(k562.apa.genes, function(x){
    tmp <- as.character(x$gene)[as.character(x$gene)%in%L]
    tmp <- length(tmp)
    return(tmp)
}))
p.values <- do.call(c,lapply(k562.apa.genes, function(x){
    tmp <- as.character(x$gene)[as.character(x$gene)%in%L]
    n1 <- length(tmp)
    m1 <- n1
    n2 <- length(as.character(x$gene)) -n1
    m2 <- length(L) - m1
    mtx <- matrix(c(n1,m1,n2,m2),c(2,2))
    res <- fisher.test(mtx)$p.value
    return(res)
}))
RBP.n <- as.data.frame(cbind(rbp.n,p.values))
RBP.n <- RBP.n[RBP.n$p.values<0.05,]
}
enrich.lineage.hepg2 <- function(L){
    rbp.n <- do.call(c,lapply(hepg2.apa.genes, function(x){
    tmp <- as.character(x$gene)[as.character(x$gene)%in%L]
    tmp <- length(tmp)
    return(tmp)
}))
p.values <- do.call(c,lapply(hepg2.apa.genes, function(x){
    tmp <- as.character(x$gene)[as.character(x$gene)%in%L]
    n1 <- length(tmp)
    m1 <- n1
    n2 <- length(as.character(x$gene)) -n1
    m2 <- length(L) - m1
    mtx <- matrix(c(n1,m1,n2,m2),c(2,2))
    res <- fisher.test(mtx)$p.value
    return(res)
}))
RBP.n <- as.data.frame(cbind(rbp.n,p.values))
RBP.n <- RBP.n[RBP.n$p.values<0.05,]
}


all.res.enriched <- lapply(all.hs.types, enrich.lineage)
all.res.mod <- lapply(all.res.enriched, function(x){
    x$RBP <- rownames(x)
    x <- x[,c('rbp.n','RBP')]
    return(x)
})

all.res <- Reduce(function(x,y)merge(x,y,by='RBP',all=T),all.res.mod)
colnames(all.res) <- c('RBP',names(all.res.mod))
rownames(all.res) <- all.res$RBP
#write.table(all.res,here('Analysis/rawdata/hs.network.raw.txt'),col.names=T,row.names=F,sep='\t',quote=F)

#hepg2============
all.res.enriched.h <- lapply(all.hs.types, enrich.lineage.hepg2)
all.res.mod.h <- lapply(all.res.enriched.h[-22], function(x){
    x$RBP <- rownames(x)
    x <- x[,c('rbp.n','RBP')]
    return(x)
})

all.res.h <- Reduce(function(x,y)merge(x,y,by='RBP',all=T),all.res.mod.h)
colnames(all.res.h) <- c('RBP',names(all.res.mod.h)[-22])
rownames(all.res.h) <- all.res.h$RBP
nodes.h <- data.frame(id=names(all.res.mod.h),n.rbps=do.call(c,lapply(all.res.mod.h, nrow)))
#save(nodes.h,all.res.h,file=here('Analysis/rawdata/hepg2.net.rawdata.RData'))
##########





#plot k562
all.res <- all.res[,-1]
all.res[is.na(all.res)] <- 0
all.res[all.res >1] <- 1
#
links <- t(all.res)
links <- links[rowSums(links)!=0,]

nodes <- data.frame(id=names(all.res.mod),n.rbps=do.call(c,lapply(all.res.mod, nrow)))

write.table(nodes,file=here('Analysis/rawdata/n.rbps.each.group.txt'),col.names=T,row.names=F,sep='\t',quote=F)
nodes <- nodes%>%mutate(size=ifelse(n.rbps>=50,10,
                                    ifelse(n.rbps>=40&n.rbps<50,9,
                                    ifelse(n.rbps >=30&n.rbps<40,8,
                                    ifelse(n.rbps>=20&n.rbps<30,7,
                                    ifelse(n.rbps>=10&n.rbps<20,6,7))))))
nodes <- nodes[nodes$n.rbps!=0,]
#
net2 <- graph_from_incidence_matrix(links)
#V(net2)$size <- c(nodes$n.rbps,rep(1,129))
#E(net2)$width <- c(nodes$n.rbps,rep(1,129))
#table(V(net2)$type)
net2.bp <- bipartite.projection(net2)
 as_incidence_matrix(net2)  %*% t(as_incidence_matrix(net2)) 

 t(as_incidence_matrix(net2)) %*%   as_incidence_matrix(net2)
V(net2.bp$proj1)$size <- nodes$n.rbps
pdf(here('test.pdf'))
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1)) 
for (layout in layouts) {
print(layout)
l <- do.call(layout, list(net2))
plot(net2, vertex.label.color=c(rep('black',nrow(links)),rep("black",ncol(links))), 
vertex.size=c(nodes$size,rep(3,ncol(links)))*3,
#vertex.size=log(nodes$n.rbps)*10, 
#edge.width =c(E(net2.bp$proj1)$weight,rep(0.5,length(E(net2.bp$proj2)$weight)))*0.2,
edge.curved=0.3,
vertex.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links)),rep('grey',length(E(net2.bp$proj2)$weight))),
vertex.label.cex=c(rep(0.5,nrow(links)),rep(0.3,ncol(links))),
vertex.frame.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links)),rep('grey',length(E(net2.bp$proj2)$weight))),
edge.arrow.mode=0,layout=l, main=layout) }
#plot(net2.bp$proj1, vertex.label.color="black", vertex.label.dist=0.1, vertex.lable=NA,vertex.size=log(nodes$n.rbps)*5,edge.curved=0.5,
#layout=layout_as_star)
dev.off()


pdf(here('Analysis/figures/network.pdf'))
plot(net2, vertex.label.color=c(rep('black',nrow(links)),rep("black",ncol(links))), 
vertex.size=c(nodes$size,rep(3,ncol(links)))*3,
#vertex.size=log(nodes$n.rbps)*10, 
#edge.width =c(E(net2.bp$proj1)$weight,rep(0.5,length(E(net2.bp$proj2)$weight)))*0.2,
edge.curved=0.2,
vertex.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links)),rep('grey',length(E(net2.bp$proj2)$weight))),
vertex.label.cex=c(rep(0.5,nrow(links)),rep(0.3,ncol(links))),
vertex.frame.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links)),rep('grey',length(E(net2.bp$proj2)$weight))),
edge.arrow.mode=0, layout=layout_with_kk)
dev.off()

#regulatory work
require(ggplot2)
require(network)
require(igraph)
require(sna)
require(ggnet)
require(ergm)
require(intergraph)
require(RColorBrewer)
require(reshape2)
te <- events.df%>%group_by(cells,detected)%>%count()
tmp <- dcast(te,cells~detected)
tmp <- all.res
tmp[is.na(tmp)] <- 0
rownames(tmp) <- tmp$RBP
tmp <- as.matrix(tmp[,-1])

g <- graph.adjacency(all.res)


bip = network(tmp,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")
ggnet2(bip, label = TRUE,edge.size = "weights",edge.alpha = 1/10)
dev.off()




d <- dist(rbp.per, method = "euclidean")
d <- as.matrix(d)

hc1 <- hclust(d, method = "complete" ,k=6)

plot(hc1, cex = 0.6, hang = -1)
pheatmap(d)
dev.off()
library(fviz)
df <- melt(d)
df <- na.omit(df)
 library("factoextra")
  library("FactoMineR")

p.z<- pheatmap(d,cluster_cols = T,show_rownames = F)

  pa.clust <- as.data.frame(cbind(d,
  cluster = cutree(p.z$tree_row,
    k = 5)))
c2 <- rownames(pa.clust[pa.clust$cluster==2,])
c1 <- rownames(pa.clust[pa.clust$cluster==1,])
c3 <- rownames(pa.clust[pa.clust$cluster==3,])
c4 <- rownames(pa.clust[pa.clust$cluster==4,])
c5 <- rownames(pa.clust[pa.clust$cluster==5,])

# Triangle heatmap to compare cohorts
ggplot(df, aes(Var1, Var2)) +
    ggtitle('Retention by cohort') +
    theme_bw() +
    xlab('Cohort') +
    ylab('Tenure (weeks)') +
    geom_tile(aes(fill = value), color='white') +
    scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab') +
    theme(axis.text.x=element_text(angle=90),
          axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'))
dev.off()
