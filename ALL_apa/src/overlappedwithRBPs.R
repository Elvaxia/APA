#overlapped number of genes when RBPs down regulated in K562/HepG2, and ALL
library(here)
library(data.table)
library(RColorBrewer)
library(igraph)
RBPs <- c('RBFOX2','NUSAP1','YBX3','MSI2','TRIP6','IGF2BP3','RCC2','SLBP','MATR3',
'PA2G4','HLTF','HNRNPA0','HNRNPAB','ILF2','KHDRBS1','NUDT21','KHSRP','DROSHA','TIA1','SFPQ','AKAP1','SUGP2','ILF3','WRN',
'BOP1','HSPD1','PHF6','GRWD1','SSRP1','EWSR1','CPSF6','NONO','SRSF1','PABPN1','CCDC124','BUD13','GEMIN5',
'DKC1','DHX30','CPSF1','UPF2','SNRNP200','HNRNPM','PCBP2','NCBP2','LARP4','TARDBP','FTO','SRPK2','QKI','FKBP4','CPEB4')

load('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP/Analysis/rawdata/all.apa.genes.2.21.RData')

genes.rbp <- k562.apa.genes['RBFOX2'][[1]]

#all shortended genes in ALL
filtered <- fread(here('/Analysis/rawdata/filtered.RUD.txt'))
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

S_genes <- unique(do.call(c,lapply(strsplit(shortened_genes,'\\:'), function(x)x[1]))) 

RBP.res<- c()
for(i in RBPs){
    RBP <- i
    RBP.res[i]<- table(as.character(k562.apa.genes[[RBP]][k562.apa.genes[[RBP]]$lable=='Lengthened',]$gene)%in%S_genes)[2]
}

RBP <-'CPSF6'
table(as.character(k562.apa.genes[[RBP]][k562.apa.genes[[RBP]]$lable=='Lengthened',]$gene)%in%S_genes)
mtx <-matrix(c(122,12,291,12),c(2,2))
fisher.test(mtx)
#p-value = 0.04148

RBP <-'TIA1'
table(as.character(k562.apa.genes[[RBP]][k562.apa.genes[[RBP]]$lable=='Lengthened',]$gene)%in%S_genes)
mtx <-matrix(c(113,5,296,5),c(2,2))
fisher.test(mtx)
#p-value = 0.15


RBP <- 'TIA1'
table(as.character(hepg2.apa.genes[[RBP]][hepg2.apa.genes[[RBP]]$lable=='Lengthened',]$gene)%in%S_genes)
mtx <-matrix(c(91,17,286,17),c(2,2))
fisher.test(mtx)
#p-value = 0.001974

RBP <- 'CPSF6'  
table(as.character(hepg2.apa.genes[[RBP]][hepg2.apa.genes[[RBP]]$lable=='Lengthened',]$gene)%in%S_genes)
mtx <-matrix(c(152,19,284,19),c(2,2))
fisher.test(mtx)
#p-value = 0.07759
library(VennDiagram)




venn.plot <- venn.diagram(
	x = list(
		CPSF6_down_K562=as.character(k562.apa.genes[['CPSF6']][k562.apa.genes[['CPSF6']]$lable=='Lengthened',]$gene),
        CPSF6_up_ALL =S_genes),
	filename = NULL,
	col = "transparent",
	fill = c("#fdae61", "#d53e4f"),
	label.col =c("#fdae61",'black','black','black',"#d53e4f",'black'),
	alpha = 0.8,
	cex = 3,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.pos = 1,
	cat.dist =  0.1,
	cat.fontfamily = "serif",
	margin = 0.1,
  height = 800, width = 800, resolution = 3000
	)
jpeg("./Analysis/figures/CPSF6.venn_jpeg.jpg",800,800);
grid.draw(venn.plot);
dev.off()

venn.plot <- venn.diagram(
	x = list(
		TIA1_down_K562=as.character(k562.apa.genes[['TIA1']][k562.apa.genes[['TIA1']]$lable=='Lengthened',]$gene),
        TIA1_up_ALL =S_genes),
	filename = NULL,
	col = "transparent",
	fill = c("#fdae61", "#d53e4f"),
	alpha = 0.5,
	cex = 3,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.pos = 3,
	cat.dist =  0.01,
	cat.fontfamily = "serif",
	margin = 0.1,
  height = 800, width = 800, resolution = 3000
	)
jpeg("./Analysis/figures/TIA1.k562.venn_jpeg.jpg",800,800);
grid.draw(venn.plot);
dev.off()


venn.plot <- venn.diagram(
	x = list(
		TIA1_down_HepG2=as.character(hepg2.apa.genes[['TIA1']][hepg2.apa.genes[['TIA1']]$lable=='Lengthened',]$gene),
        TIA1_up_ALL =S_genes),
	filename = NULL,
	col = "transparent",
	fill = c("#fdae61", "#d53e4f"),
	alpha = 0.50,,
	cex = 3,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.pos = 3,
	cat.dist =  0.01,
	cat.fontfamily = "serif",
	margin = 0.1,
  height = 800, width = 800, resolution = 3000
	)
jpeg("./Analysis/figures/TIA1.venn_jpeg.jpg",800,800);
grid.draw(venn.plot);
dev.off()

venn.plot <- venn.diagram(
	x = list(
		CPSF6_down_HepG2=as.character(hepg2.apa.genes[['CPSF6']][hepg2.apa.genes[['CPSF6']]$lable=='Lengthened',]$gene),
        CPSF6_up_ALL =S_genes),
	filename = NULL,
	col = "transparent",
	fill = c("#fdae61", "#d53e4f"),
	alpha = 0.50,,
	cex = 3,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.pos = 3,
	cat.dist =  0.01,
	cat.fontfamily = "serif",
	margin = 0.1,
  height = 800, width = 800, resolution = 3000
	)
jpeg("./Analysis/figures/CPSF6.hepg2.venn_jpeg.jpg",800,800);
grid.draw(venn.plot);
dev.off()

##%genes terget by RBPs==================================================================
library(dplyr)
load('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP/Analysis/rawdata/all.apa.genes.2.21.RData')
de.genes <- c(L_genes,S_genes)
events.df <- data.frame(gene=de.genes,type=c(rep('Lengthened',length(L_genes)),rep('Shortened',length(S_genes))))

all.apa.genes.k562 <- unique(do.call(c,lapply(k562.apa.genes, function(x)as.character(x$gene))))
all.apa.genes.hepg2 <- unique(do.call(c,lapply(hepg2.apa.genes, function(x)as.character(x$gene))))

length(de.genes[de.genes%in%all.apa.genes.k562])/length(de.genes)
length(de.genes[de.genes%in%all.apa.genes.hepg2])/length(de.genes)

##network==================================================================

enrich.ALL <- function(L,cellline){
    rbp.n <- do.call(c,lapply(cellline, function(x){
    tmp <- as.character(x[x$lable=='Lengthened',]$gene)[as.character(x[x$lable=='Lengthened',]$gene)%in%L]	
    tmp <- length(tmp)
    return(tmp)
}))
p.values <- do.call(c,lapply(cellline, function(x){
    tmp <- as.character(x[x$lable=='Lengthened',]$gene)[as.character(x[x$lable=='Lengthened',]$gene)%in%L]
    n1 <- length(tmp)
    m1 <- n1
    n2 <- length(as.character(x[x$lable=='Lengthened',]$gene)) -n1
    m2 <- length(L) - m1
    mtx <- matrix(c(n1,m1,n2,m2),c(2,2))
    res <- fisher.test(mtx)$p.value
    return(res)
}))
RBP.n <- as.data.frame(cbind(rbp.n,p.values))
RBP.n <- RBP.n[RBP.n$p.values<0.05,]
}


rbps.ALL.k562 <- enrich.ALL(S_genes,k562.apa.genes)


#use DE rbps
de.rbps.all <- fread(here('Analysis/rawdata/RBPs.expression.txt'))

use.rbps <- intersect(rownames(rbps.ALL.k562),de.rbps.all$gene_name)
rbps.ALL.k562 <- rbps.ALL.k562[use.rbps,]
hs.network <- fread('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_hs/Analysis/rawdata/hs.network.raw.txt')


rbps.ALL.k562 <- merge(rbps.ALL.k562,hs.network,by.x='row.names',by.y='RBP',all=T)

rownames(rbps.ALL.k562) <- rbps.ALL.k562$Row.names
rbps.ALL.k562 <- rbps.ALL.k562[,-c(1:2)]

colnames(rbps.ALL.k562)[1] <- 'ALL'

rbps.ALL.k562[is.na(rbps.ALL.k562)] <- 0
rbps.ALL.k562[rbps.ALL.k562 >0] <- 1

links <- t(rbps.ALL.k562)
links <- links[rowSums(links)!=0,]

nodes.hs <- fread('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_hs/Analysis/rawdata/n.rbps.each.group.txt')


nodes <- data.frame(id=c('ALL',nodes.hs$id),n.rbps=c(length(use.rbps),nodes.hs$n.rbps))
nodes <- nodes%>%mutate(size=ifelse(n.rbps>=50,10,
                                    ifelse(n.rbps>=40&n.rbps<50,9,
                                    ifelse(n.rbps >=30&n.rbps<40,8,
                                    ifelse(n.rbps>=20&n.rbps<30,7,
                                    ifelse(n.rbps>=10&n.rbps<20,6,7))))))
nodes <- nodes[nodes$n.rbps!=0,]
#only use lymphoid cell types
links.l <- links[c('ALL','B','Lymphoid','CLP','MPP','HSC'),]
links.l <- links.l[,colSums(links.l)!=0]
nodes.l <- nodes[nodes$id%in%c('ALL','B','Lymphoid','CLP','MPP','HSC'),]

#
net2 <- graph_from_incidence_matrix(links.l)
#V(net2)$size <- c(nodes$n.rbps,rep(1,129))
#E(net2)$width <- c(nodes$n.rbps,rep(1,129))
#table(V(net2)$type)
net2.bp <- bipartite.projection(net2)
 as_incidence_matrix(net2)  %*% t(as_incidence_matrix(net2)) 

 t(as_incidence_matrix(net2)) %*%   as_incidence_matrix(net2)

V(net2.bp$proj1)$size <- nodes.l$n.rbps

pdf(here('test.pdf'))
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1)) 
for (layout in layouts) {
print(layout)
l <- do.call(layout, list(net2))
plot(net2, vertex.label.color='black', 
vertex.size=c(nodes.l$size,rep(3,ncol(links.l)))*3,
#vertex.size=log(nodes$n.rbps)*10, 
#edge.width =c(E(net2.bp$proj1)$weight,rep(0.5,length(E(net2.bp$proj2)$weight)))*0.2,
edge.curved=0.3,
vertex.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
vertex.label.cex=c(rep(0.5,nrow(links.l)),rep(0.3,ncol(links.l))),
vertex.frame.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
edge.arrow.mode=0,layout=l, main=layout) }
#plot(net2.bp$proj1, vertex.label.color="black", vertex.label.dist=0.1, vertex.lable=NA,vertex.size=log(nodes$n.rbps)*5,edge.curved=0.5,
#layout=layout_as_star)
dev.off()





pdf(here('Analysis/figures/netrowk.all.lym.pdf'))
plot(net2, vertex.label.color='black', 
vertex.size=c(nodes.l$size*1.2,rep(3.8,ncol(links.l)))*3.7,
#vertex.size=log(nodes$n.rbps)*10, 
#edge.width =c(E(net2.bp$proj1)$weight,rep(0.5,length(E(net2.bp$proj2)$weight)))*0.2,
edge.curved=0.3,
vertex.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
vertex.label.cex=c(rep(0.8,nrow(links.l)),rep(0.45,ncol(links.l))),
vertex.frame.color=c(colorRampPalette(brewer.pal(10,"Set3")[-9])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
edge.arrow.mode=0,layout=layout_with_kk)
dev.off()




#==========HepG2============================================
rbps.ALL.hepg2 <- enrich.ALL(S_genes,hepg2.apa.genes)
use.rbps.h <- intersect(rownames(rbps.ALL.hepg2),de.rbps.all$gene_name)
rbps.ALL.hepg2 <- rbps.ALL.hepg2[use.rbps.h,]

load('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_hs/Analysis/rawdata/hepg2.net.rawdata.RData')
rbps.ALL.hepg2 <- merge(rbps.ALL.hepg2,all.res.h,by.x='row.names',by.y='RBP',all=T)

rownames(rbps.ALL.hepg2) <- rbps.ALL.hepg2$Row.names
rbps.ALL.hepg2 <- rbps.ALL.hepg2[,-c(1:2)]

colnames(rbps.ALL.hepg2)[1] <- 'ALL'

rbps.ALL.hepg2[is.na(rbps.ALL.hepg2)] <- 0
rbps.ALL.hepg2[rbps.ALL.hepg2 >0] <- 1

links <- t(rbps.ALL.hepg2)
links <- links[rowSums(links)!=0,]

nodes <- data.frame(id=c('ALL',as.character(nodes.h$id)),n.rbps=c(length(use.rbps.h),nodes.h$n.rbps))
nodes <- nodes%>%mutate(size=ifelse(n.rbps>=50,10,
                                    ifelse(n.rbps>=40&n.rbps<50,9,
                                    ifelse(n.rbps >=30&n.rbps<40,8,
                                    ifelse(n.rbps>=20&n.rbps<30,7,
                                    ifelse(n.rbps>=10&n.rbps<20,6,7))))))
nodes <- nodes[nodes$n.rbps!=0,]

#only use lymphoid cell types
links.l <- links[c('ALL','Lymphoid','CLP','MPP','HSC'),]
links.l <- links.l[,colSums(links.l)!=0]
nodes.l <- nodes[nodes$id%in%c('ALL','Lymphoid','CLP','MPP','HSC'),]

net2 <- graph_from_incidence_matrix(links.l)
#V(net2)$size <- c(nodes$n.rbps,rep(1,129))
#E(net2)$width <- c(nodes$n.rbps,rep(1,129))
#table(V(net2)$type)
net2.bp <- bipartite.projection(net2)
 as_incidence_matrix(net2)  %*% t(as_incidence_matrix(net2)) 

 t(as_incidence_matrix(net2)) %*%   as_incidence_matrix(net2)

V(net2.bp$proj1)$size <- nodes.l$n.rbps


pdf(here('Analysis/figures/netrowk.all.lym.hepg2.pdf'))
plot(net2, vertex.label.color='black', 
vertex.size=c(nodes.l$size*1.2,rep(3.8,ncol(links.l)))*3.7,
#vertex.size=log(nodes$n.rbps)*10, 
#edge.width =c(E(net2.bp$proj1)$weight,rep(0.5,length(E(net2.bp$proj2)$weight)))*0.2,
edge.curved=0.3,
vertex.color=c(colorRampPalette(brewer.pal(10,"Set3")[c(-2,-9)])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
vertex.label.cex=c(rep(0.8,nrow(links.l)),rep(0.45,ncol(links.l))),
vertex.frame.color=c(colorRampPalette(brewer.pal(10,"Set3")[c(-2,-9)])(nrow(links.l)),rep('grey',length(E(net2.bp$proj2)))),
edge.arrow.mode=0,layout=layout_with_kk)
dev.off()
