#detected in many(>5celltype),detetcted in all celltye(17),detetcted in some(1<x<=5), detetcted in single celltype
#rbp apa genes
library(here)
library(circlize)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
load('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP/Analysis/rawdata/all.apa.genes.2.21.RData')
load(here('Analysis/rawdata/hs.group.detected.RData'))

te <- events.df%>%group_by(cells,detected)%>%count()
grid.col <- c(Erythroid = '#beaed4', Lymhoid='#7fc97f',Myeloid='#fdc086',Progenitors='#ffff99',
    Single='#386cb0',HSC = '#8dd3c7',CMP='#ffffb3',GMP='#ffffb3',MEP='#ffffb3',
    CLP='#fccde5',MPP='#fccde5',EB='#bebada',MK='#bebada',DC='#fb8072',
    Mp='#fb8072',Mo='#80b1d3',Neu='#80b1d3',CD4='#fdb462',CD8='#fdb462',B='#33a02c',NK='#33a02c')
    
pdf(here('Analysis/figures/circular.group.pdf'),8,8)
chordDiagram(te,order=c('HSC','Mo','EB','NK','GMP','MK','DC','Mp','CMP',
    'Neu','MEP','CD4','CD8','CLP','B','MPP','Lymphoid','Single',
    'Progenitors','Myeloid','Erythoid'),grid.col = grid.col)
dev.off()

#hsc genes covered by RBP
all.apa.genes.k562 <- unique(do.call(c,lapply(k562.apa.genes, function(x)as.character(x$gene))))

all.apa.genes.hepg2 <- unique(do.call(c,lapply(hepg2.apa.genes, function(x)as.character(x$gene))))
events.df <- events.df%>%mutate(coveredBYrbp = ifelse(gene %in%all.apa.genes.k562,'Covered','Uncovered'))
events.df <- events.df%>%mutate(coveredBYrbp.h = ifelse(gene %in%all.apa.genes.hepg2,'Covered.h','Uncovered.h'))

use.df.k562 <- events.df %>% group_by(cells,coveredBYrbp)%>%count()%>%group_by(cells)%>%
            mutate(pro=prop.table(n)*100) %>% ungroup%>%filter(coveredBYrbp=='Covered')


use.df.hepg2 <- events.df %>% group_by(cells,coveredBYrbp.h)%>%count()%>%group_by(cells)%>%
            mutate(pro=prop.table(n)*100) %>% ungroup%>%filter(coveredBYrbp.h=='Covered.h')

use.df <- data.frame(cells=use.df.k562$cells,K562=use.df.k562$pro,HepG2=use.df.hepg2$pro)
#use groups

use.gr.k562 <- events.df %>% group_by(detected,coveredBYrbp)%>%count()%>%group_by(detected)%>%
            mutate(pro=prop.table(n)*100) %>% ungroup%>%filter(coveredBYrbp=='Covered')
use.gr.hepg2 <- events.df %>% group_by(detected,coveredBYrbp.h)%>%count()%>%group_by(detected)%>%
            mutate(pro=prop.table(n)*100) %>% ungroup%>%filter(coveredBYrbp.h=='Covered.h')

use.gr <- data.frame(Group = use.gr.k562$detected,K562=use.gr.k562$pro,HepG2=use.gr.hepg2$pro)

names <- c('HSC','GMP','CMP','MEP','CLP','MPP','EB','MK','UVEC','DC','Mp','Mo','Neu','CD4','CD8','B','NK')


pdf(here('Analysis/figures/per.gene.targeted.byRBP.pdf'),height=2.5,width=8)
ggplot(melt(use.df),aes(x=cells,y=variable))+geom_tile(aes(fill = value))+
scale_fill_continuous(low='#fde0dd',high='#c51b8a')+
scale_x_discrete(limits=names)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        #legend.title=element_blank(),
                        legend.position='top',
                        legend.text=element_text(size=12),
                        axis.title = element_text(size=22),
                        axis.text.x = element_text(angle = -30),
                        axis.text = element_text(size=12))+
labs(x='',y='',fill='% Genes targeted by RBP') 
dev.off()

pdf(here('Analysis/figures/per.gene.group.targeted.byRBP.pdf'),height=2.5,width=6)
ggplot(melt(use.gr),aes(x=Group,y=variable))+geom_tile(aes(fill = value))+
scale_fill_continuous(low='#fde0dd',high='#c51b8a')+
#scale_x_discrete(limits=names)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        #legend.title=element_blank(),
                        legend.position='top',
                        legend.text=element_text(size=12),
                        axis.title = element_text(size=22),
                        axis.text.x = element_text(angle = -30),
                        axis.text = element_text(size=12))+
labs(x='',y='',fill='% Genes targeted by RBP') 
dev.off()


#RBP expression in celltypes
#group1: detected 20% genes in single detected group

#group2: detetcted in all celltypes(17 celltypes)
#group3: detetced in some(1,5)
#group4: detetced in single celltype



events.df <- events.df%>% mutate(gene = gsub(':.*','',pa))
tmp <- intersect(k562.apa.genes[[1]]$gene, events.df$gene)
x <- 2
groups <- do.call(rbind,lapply(k562.apa.genes, function(x){
    n.groups <- table(events.df[events.df$gene%in%x$gene,]$detected)
    n.genes <- length(x$gene[x$gene%in%unique(events.df$gene)]) / length(x$gene)
    df <- data.frame(ncells=sum(n.groups!=0),ngenes=n.genes)
    return(df)
}))

use.groups <- groups%>%muate(ifelse(ngenes>))

#cluster of RBPs 
require(ggplot2)
require(network)
require(igraph)
require(sna)
require(ggnet)
require(ergm)
require(intergraph)
require(RColorBrewer)
require(reshape2)

tmp <- dcast(te,cells~detected)
tmp[is.na(tmp)] <- 0
rownames(tmp) <- tmp$cells
tmp <- as.matrix(tmp[,-1])


                 row.names = letters[1:4])
bip = network(tmp,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")
ggnet2(bip, label = TRUE,edge.size = "weights",edge.alpha = 1/10)
dev.off()

#cluster RBPs might contribute to HSC differetiation
#number of differential APA
library(rgexf)
library(RCytoscape)