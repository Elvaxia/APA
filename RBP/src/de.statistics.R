#plot filter statistics=====================================================
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP')
library(here)
library(data.table)
library(dplyr)
library(ggplot2)
require("ggrepel")
load(here('Analysis/rawdata/rbp.testres.2.21.RData'))


stat.df <- function(rbp.DE){
  n.apa.events <- do.call(c,lapply(rbp.DE,nrow))
n.de.apa <- do.call(c,lapply(rbp.DE,function(x){
  nrow(x[p_value < 0.05& abs(log2FC) >= 0.5,])
}))
n.shortened <- do.call(c,lapply(rbp.DE,function(x){
    nrow(x[ p_value <0.05 & log2FC<= -0.5,])
}))
n.lengthened <- do.call(c,lapply(rbp.DE,function(x){
    nrow(x[p_value <0.05 & log2FC >= 0.5,])
}))
df <- as.data.frame(cbind(n.apa.events,n.de.apa,n.shortened,n.lengthened))
df$Names <- rownames(df)
return(df)
}



hepg2.stat <- stat.df(hepg2.de)
k562.stat <- stat.df(new.k562.de)

##
apa.factors <- c('CPSF6','CPSF7','CSTF2T','CSTF2','CELF1','CPSF1','CPSF2','CPSF3','CPSF4','FIP1L1',
'CSTF1','CSTF3','NUDT21','PCF11','CLP1','CPSF3L','SYMPK','WDR33','RBBP6','CPEB1','PAPOLA',
'PAPOLG','PAPD4','PAPD5','PAPD7','PABPC4','PABPN1','PABPC')
regulators <- fread(here('regulators_apa_nature_2018.txt'))
##6 APA factors in RBPs RNA-seq


pdf(here('Analysis/figures/filter.statistic.hepg2.pdf'),8,8)
df <- hepg2.stat
df <- df %>% mutate(Colors = ifelse(Names%in%apa.factors,'APA factors','Others'))
ggplot(df,aes(x=n.apa.events,y=n.de.apa,label=Names))+geom_point(aes(color = factor(Colors)),size=5)+
geom_label_repel(aes(label=ifelse(Colors == 'APA factors' ,as.character(Names),'')),
color ='red',size=5)+
geom_label_repel(aes(label=ifelse(n.de.apa>=480  ,as.character(Names),'')), color = 'grey',
size=5) +ggtitle('HepG2')+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
                  legend.position = 'top', legend.title = element_blank(),
                  legend.text=element_text(size=16),
                  axis.title.x = element_text(size=25),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=25),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
labs(x= 'Number of total APA events', y='Number of differential APA events') + xlim(0,9000)
dev.off()


pdf(here('Analysis/figures/filter.statistic.k562.pdf'),8,8)
df <- k562.stat
df <- df %>% mutate(Colors = ifelse(Names%in%apa.factors,'APA factors','Others'))
ggplot(df,aes(x=n.apa.events,y=n.de.apa,label=Names))+geom_point(aes(color = factor(Colors)),size=5)+
geom_label_repel(aes(label=ifelse(n.de.apa>=480&Colors!='APA factors'  ,as.character(Names),'')), color = 'grey',vjust = 1,
size=5,hjust = -0.5)+
geom_label_repel(aes(label=ifelse(Colors == 'APA factors' ,as.character(Names),'')),
color ='red',size=5,vjust= 0.5)+
 ggtitle('K562')+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
                  legend.position = 'top', legend.title = element_blank(),
                  legend.text=element_text(size=16),
                  axis.title.x = element_text(size=25),
                  axis.text.x =element_text(size=16),
                  axis.title.y = element_text(size=25),
                  axis.text.y  = element_text(size=16),
                  panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))+
labs(x= 'Number of total APA events', y='Number of differential APA events') + xlim(0,9000)
dev.off()



s.l.col <- c("#fc8d62","#66c2a5")
names(s.l.col) <- c("Lengthened", "Shortened")


df <- k562.stat
colnames(df)[3:4] <- c("Shortened","Lengthened")
df$Percent <- round(df$Shortened/df$n.de.apa,2) *100
name.ordered <- df[order(df$n.de.apa,decreasing=F),]$Names
df <- df%>% mutate(Type = ifelse(Percent >50,"Shortened preference","Lengthened preference"))
data.use <- melt(df[,3:5],id='Names')
data.use <- merge(data.use,df[,c(5,7)],by='Names')
data.use$cell <-'K562'

df2 <- hepg2.stat
colnames(df2)[3:4] <- c("Shortened","Lengthened")
df2$Percent <- round(df2$Shortened/df2$n.de.apa,2) *100
df2 <- df2%>% mutate(Type = ifelse(Percent >50,"Shortened preference","Lengthened preference"))
data.use2 <- melt(df2[,3:5],id='Names')
data.use2 <- merge(data.use2,df2[,c(5,7)],by='Names')
data.use2$cell <- 'HepG2'

data.use.all <- rbind(data.use,data.use2)

both.l <- intersect(df[df$Type=='Lengthened preference',]$Names,df2[df2$Type=='Lengthened preference',]$Names)
both.s <- intersect(df[df$Type=='Shortened preference',]$Names,df2[df2$Type=='Shortened preference',]$Names)
consistance_RBPs <- c(both.l,both.s)


name.ordered <- df[order(df$n.de.apa,decreasing=F),]$Names
name.ordered.l <- intersect(name.ordered,both.l)
name.ordered.s <- intersect(name.ordered,both.s)

pdf(here::here('Analysis/figures/k562.hepg2.l.pdf'),width=4,height=8)
ggplot(data.use.all[data.use.all$Names%in%consistance_RBPs,],aes(x=Names,y=value,fill=variable)) + 
geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
facet_grid(. ~cell,scales='free',space='free') +
scale_x_discrete(limits=name.ordered.l)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown',title='Lengthened preference') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
dev.off()
pdf(here::here('Analysis/figures/k562.hepg2.s.pdf'),width=4,height=5)
ggplot(data.use.all[data.use.all$Names%in%consistance_RBPs,],aes(x=Names,y=value,fill=variable)) + 
geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
facet_grid(. ~cell,scales='free',space='free') +
scale_x_discrete(limits=name.ordered.s)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown',title='Shortened preference') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
dev.off()
pdf(here::here('Analysis/figures/k562.hepg2.incons.pdf'),width=4,height=12)
ggplot(data.use.all[!data.use.all$Names%in%consistance_RBPs,],aes(x=Names,y=value,fill=variable)) + 
geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
facet_grid(. ~cell,scales='free',space='free') +
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown',title='Inconsistent') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
dev.off()

pdf(here('Analysis/figures/s.l.hist.statistic.k562.pdf'),height=8,width=4)
df <- k562.stat
colnames(df)[3:4] <- c("Shortened","Lengthened")
data.use <- melt(df[,3:5],id='Names')
name.ordered <- df[order(df$n.de.apa,decreasing=F),]$Names
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[1:50])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()


ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[51:100])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[101:150])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[151:211])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
dev.off()



pdf(here::here('Analysis/figures/s.l.hist.statistic.hepg2.pdf'),height=8,width=4)
df <- hepg2.stat
colnames(df)[3:4] <- c("Shortened","Lengthened")
name.ordered <- df[order(df$n.de.apa,decreasing=F),]$Names
data.use <- melt(df[,3:5],id='Names')
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[1:50])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[51:100])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[101:150])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
ggplot(data.use,aes(x=Names,y=value,fill=variable)) + geom_bar(stat="identity") +
scale_fill_manual(values=s.l.col)+
scale_x_discrete(limits=name.ordered[151:199])+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=16),
                        axis.text = element_text(size=8))+
labs(x='',y='Number of genes after each RBP knockdown') +
scale_y_continuous(position = "right",expand = c(0, 0),limits=c(0,700)) +
coord_flip()
dev.off()


library(VennDiagram)
#overlaped RBPs in HepG2 and K562

#overlapped RBPs both target more genes shortend or lengthend in HepG2 and K562 
hepg2.stat <- hepg2.stat%>% mutate(percentageOFs = n.shortened/n.de.apa)
k562.stat <- k562.stat%>% mutate(percentageOFs = n.shortened/n.de.apa)

hepg2.all <- hepg2.stat$Names
k562.all <- k562.stat$Names
length(hepg2.all)
# 199
length(k562.all)
#211
length(unique(c(hepg2.all,k562.all)))
#239
both.s <- intersect(k562.stat[k562.stat$percentageOFs >= 0.5,]$Names ,hepg2.stat[hepg2.stat$percentageOFs > 0.5,]$Names)
both.l <- intersect(k562.stat[k562.stat$percentageOFs < 0.5,]$Names ,hepg2.stat[hepg2.stat$percentageOFs < 0.5,]$Names)
length(both.s)
#31
length(both.l)
#75
length(K562.lengthened)
#127
length(HepG2.lengthened)
#130
length(HepG2.shortened)
#69
length(K562.shortened)
#84

venn.plot <- venn.diagram(
	x = list(
		K562.lengthened = k562.stat[k562.stat$percentageOFs < 0.5,]$Names,
    HepG2.lengthened = hepg2.stat[hepg2.stat$percentageOFs < 0.5,]$Names,
		K562.shortened = k562.stat[k562.stat$percentageOFs >= 0.5,]$Names,
    HepG2.shortened = hepg2.stat[hepg2.stat$percentageOFs >= 0.5,]$Names
		),
	filename = NULL,
	col = "transparent",
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
	alpha = 0.50,
	label.col = c("orange", "white", "darkorchid4", "white", 
	"white", "white", "white", "white", "darkblue", "white", 
	"white", "white", "white", "darkgreen", "white"),
	cex = 2,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
	cat.cex = 2,
	cat.pos = 0,
	cat.dist =  0.1,
	cat.fontfamily = "serif",
	margin = 0.1,
  height = 800, width = 800, resolution = 600
	)
jpeg("./Analysis/figures/venn_jpeg.jpg",800,800);
grid.draw(venn.plot);
dev.off()


#==========================================================================================================
rbp <- 'CPSF6'
library(ggplot2)

pdf(here('Analysis/figures/cpsf6.k562.significant.statistics.pdf'),height=6,width=8)
rbp.DE.dt <- new.k562.de[[rbp]]
rbp.DE.dt[rbp.DE.dt$kd.rud.mean >1,]$kd.rud.mean <- 1
rbp.DE.dt[rbp.DE.dt$ctrl.rud.mean >1,]$ctrl.rud.mean <- 1
rbp.DE.dt <- na.omit(rbp.DE.dt)
rbp.DE.dt <- rbp.DE.dt %>% mutate(res = ifelse(p_value <0.05,
                                                ifelse(log2FC >= 0.5,'Lengthened',
                                                ifelse(log2FC <= -0.5,'Shortened','Not significant')),'Not significant'))
ggplot(rbp.DE.dt,aes(x= log2FC,y= -log10(p_value),fill=factor(res))) + 
  geom_point(aes(col=factor(res))) +
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

ggplot(rbp.DE.dt[rbp.DE.dt$res!='Not significant',],aes(x=ctrl.rud.mean,y=kd.rud.mean))+
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
                        xlim(c(0,1))+ylim(c(0,1))+labs(x= 'RUCP of control', y='RUCP of knockdown')+geom_abline(slope=1, intercept=0) +
 annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(rbp.DE.dt$res)[1]),color='red',size=6) + 
  annotate("text", x = 0.8, y=0.1, label = paste0('Shortened',' ',table(rbp.DE.dt$res)[3]),color ='blue',size=6)
dev.off()
library(hexbin)
pdf(here('Analysis/figures/cpsf6.k562.hexbin.statistics.pdf'),8,6)
ggplot(rbp.DE.dt[rbp.DE.dt$res!='Not significant',],aes(x=ctrl.rud.mean,y=kd.rud.mean)) +
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
  labs(x= 'RUCP of control', y='RUCP of knockdown')+geom_abline(slope=1, intercept=0) +
  annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(rbp.DE.dt$res)[1]),color='red',size=5.5) + 
  annotate("text", x = 0.8, y=0.1, label = paste0('Shortened',' ',table(rbp.DE.dt$res)[3]),color ='blue',size=5.5)
dev.off()


pdf(here('Analysis/figures/cpsf6.hepg2.significant.statistics.pdf'),height=6,width=8)
rbp.DE.dt <- hepg2.de[[rbp]]
rbp.DE.dt[rbp.DE.dt$kd.rud.mean >1,]$kd.rud.mean <- 1
rbp.DE.dt[rbp.DE.dt$ctrl.rud.mean >1,]$ctrl.rud.mean <- 1
rbp.DE.dt <- na.omit(rbp.DE.dt)
rbp.DE.dt <- rbp.DE.dt %>% mutate(res = ifelse(p_value <0.05,
                                                ifelse(log2FC >= 0.5,'Lengthened',
                                                ifelse(log2FC <= -0.5,'Shortened','Not significant')),'Not significant'))
ggplot(rbp.DE.dt,aes(x= log2FC,y= -log10(p_value),fill=factor(res))) + geom_point(aes(col=factor(res))) +
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
ggplot(rbp.DE.dt[rbp.DE.dt$res!='Not significant',],aes(x=ctrl.rud.mean,y=kd.rud.mean))+geom_point(size=6)+
theme(panel.background=element_blank(),
                        panel.grid=element_blank(),
                        axis.line=element_line(),
                        plot.margin = margin(t=0, b=0.5, r=2, l=0.5, "cm"),
                        legend.title=element_blank(),
                        legend.position="bottom",
                        legend.text=element_text(size=8),
                        axis.title = element_text(size=22),
                        axis.text = element_text(size=12))+
                        xlim(c(0,1))+ylim(c(0,1))+labs(x= 'RUCP of control', y='RUCP of knockdown')+geom_abline(slope=1, intercept=0) +
 annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(rbp.DE.dt$res)[1]),color='red',size=6) + 
  annotate("text", x = 0.8, y=0.1, label = paste0('Shortened',' ',table(rbp.DE.dt$res)[3]),color ='blue',size=6)
dev.off()


pdf(here('Analysis/figures/cpsf6.k562.hexbin.hep.statistics.pdf'),8,6)
ggplot(rbp.DE.dt[rbp.DE.dt$res!='Not significant',],aes(x=ctrl.rud.mean,y=kd.rud.mean)) +
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
  labs(x= 'RUCP of control', y='RUCP of knockdown')+geom_abline(slope=1, intercept=0) +
  annotate("text", x=0.2, y=0.8, label= paste0('Lengthened',' ',table(rbp.DE.dt$res)[1]),color='red',size=5.5) + 
  annotate("text", x = 0.8, y=0.1, label = paste0('Shortened',' ',table(rbp.DE.dt$res)[3]),color ='blue',size=5.5)
dev.off()
#===================='CPSF6'===================================================
library(ggpubr)
CPSF6.k562 <- new.k562.de[['CPSF6']]
CPSF6.hepg2 <- hepg2.de[['CPSF6']]
data <- merge(CPSF6.k562[,c('delta_RUD','id')],CPSF6.hepg2[,c('delta_RUD','id')],by='id',all=F)
colnames(data) <- c('ID','K562','HepG2')


QuadStac <- function(data, x, y, P) {
  setDT(data)
  qu1 <- data[data[[x]] > 0 & data[[y]] > 0, ]
  qu1.x <- mean(range(qu1[[x]]))
  qu1.y <- mean(range(qu1[[y]]))
  qu1.lab <- round(nrow(qu1)/nrow(data)*100, 2)
  res1 <- c(qu1.x, qu1.y, qu1.lab)
  
  qu1 <- data[data[[x]] < 0 & data[[y]] > 0, ]
  qu1.x <- mean(range(qu1[[x]]))
  qu1.y <- mean(range(qu1[[y]]))
  qu1.lab <- round(nrow(qu1)/nrow(data)*100, 2)
  res2 <- c(qu1.x, qu1.y, qu1.lab)
  
  qu1 <- data[data[[x]] < 0 & data[[y]] < 0, ]
  qu1.x <- mean(range(qu1[[x]]))
  qu1.y <- mean(range(qu1[[y]]))
  qu1.lab <- round(nrow(qu1)/nrow(data)*100, 2)
  res3 <- c(qu1.x, qu1.y, qu1.lab)
  
  qu1 <- data[data[[x]] > 0 & data[[y]] < 0, ]
  qu1.x <- mean(range(qu1[[x]]))
  qu1.y <- mean(range(qu1[[y]]))
  qu1.lab <- round(nrow(qu1)/nrow(data)*100, 2)
  res4 <- c(qu1.x, qu1.y, qu1.lab)
  
  Lab <- list(res1, res2, res3, res4)
  
 colnames(data)[which(colnames(data) == x)] <- "xc"
  colnames(data)[which(colnames(data) == y)] <- "yc"
  
  pdf(here('Analysis/figures/cpsf6.cor.hepg2.k562.pdf'),8,6)
  ggplot(data = data, aes(x = xc, y = yc)) + 
    geom_point(size = .1)+
    stat_density2d(geom = "raster", aes(fill = ..density.., alpha = ..density..), contour = FALSE) +
    scale_fill_viridis_c(guide = FALSE) +
    scale_alpha_continuous(guide = "none", range = c(0, 1))+
    labs(x = x, y = y, title = P, caption = "Delta RUCP", subtitle = paste("N =", nrow(data)))+
    annotate("text", col = "red", size = 6, 
             x = c(Lab[[1]][1], Lab[[2]][1], Lab[[3]][1], Lab[[4]][1]), 
             y = c(Lab[[1]][2], Lab[[2]][2], Lab[[3]][2], Lab[[4]][2]), 
             label = c(Lab[[1]][3], Lab[[2]][3], Lab[[3]][3], Lab[[4]][3])) + 
    theme_classic() + 
    geom_hline(yintercept = 0, lty = 2) + 
    geom_vline(xintercept = 0, lty = 2) + 
    geom_smooth(method = "lm", lwd = 0.5) + 
stat_cor(data = data[, .(xc, yc)], method = "pearson", 
             size = 6, col = "red", 
             label.y.npc = "top", label.x.npc = "center") + 
    theme(panel.background = element_blank(), 
          panel.grid = element_line(colour = "grey90"), 
          axis.title = element_text(size = 22), 
          axis.text = element_text(size = 16), legend.position = "none")
dev.off()
}



QuadStac(data = data, x = "K562", y = "HepG2", P = 'CPSF6')
