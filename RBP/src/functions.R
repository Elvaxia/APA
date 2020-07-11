#function for ALL
#=====================================================================
###sashimiplot
##=====================================================================
library(GenomicFeatures)
library(grid)
library(RColorBrewer)
library(Gviz)
library(data.table)
library(here)
options(ucscChromosomeNames=FALSE)
#load gtf
gtrack <- GenomeAxisTrack()
txdb <- loadDb('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/GRCh38.93.txdb')
k562<- fread(here('RBP.sample.info.txt'),header=T)
hepg2 <- fread(here('RBP.sample.HepG2.info.txt'),header=T)
sample.info <- hepg2
sample.info <- k562
#set colors for celltypes

my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
names(my_colors) <- c('KD','Ctrl')


fig.path = here('Analysis/figures/sashimiplot/HepG2')
fig.path = here('Analysis/figures/sashimiplot/K562')
#plot for group 
sashimi.plot <- function(id,rbp){
  sample.info <- sample.info[RBP==rbp,]
  sample.info <- sample.info[order(sample.info$Group),]
  #prepare bam files
  align.list <- list()
  for(s in sample.info$sample){
    s.name <- s
    s.bam <- sample.info[sample==s,]$Bam
    color <- sample.info[sample==s,]$Group
    s.align <- AlignmentsTrack(s.bam,
                              isPaired = TRUE,
                              name=s.name,
                              fill.coverage = my_colors[color] )
    align.list[[s.name]] <- s.align
  }

  region <- strsplit(id,':')[[1]]
  from = as.numeric(region[4])
  to = as.numeric(region[6])
  chrs <- region[2]
  aligns <- align.list
  n <- length(aligns)
  if(!dir.exists(file.path(fig.path,rbp))){
    dir.create(file.path(fig.path,rbp))}
  grtrack <- GeneRegionTrack(txdb, 
                             name="Gene",
                             chromosome=chrs,
                             from = from ,
                             to = to )
  ht <- HighlightTrack(trackList=gtrack, 
                       start=as.numeric(region[5]), 
                       width=2, 
                       chromosome=chrs)
  pdf(file.path(fig.path,rbp,paste0(id,'2.29.pdf')),width = 8, height = 8 )
  plotTracks(c(ht,grtrack,aligns),
               from = from ,
               to = to ,
               chromosome = chrs,
               sizes = c(0.4,0.6,rep(1,n)),
               #when error :Too many stacks to draw. ,
               #remenber expand grtrack size.
               extend.right = 50,
               extend.left = 50,
               type="coverage",
               showSampleNames=TRUE,
               main = id,
               cex.main = 1,
               labelPos="above")
    dev.off()
}

sashimi.plot.log <- function(id,rbp){
  sample.info <- sample.info[RBP==rbp,]
  sample.info <- sample.info[order(sample.info$Group),]
  #prepare bam files
  align.list <- list()
  for(s in sample.info$sample){
    s.name <- s
    s.bam <- sample.info[sample==s,]$Bam
    color <- sample.info[sample==s,]$Group
    s.align <- AlignmentsTrack(s.bam,
                              isPaired = TRUE,
                              name=s.name,
                              fill.coverage = my_colors[color] )
    align.list[[s.name]] <- s.align
  }

  region <- strsplit(id,':')[[1]]
  from = as.numeric(region[4])
  to = as.numeric(region[6])
  chrs <- region[2]
  aligns <- align.list
  n <- length(aligns)
  if(!dir.exists(file.path(fig.path,rbp))){
    dir.create(file.path(fig.path,rbp))}
  grtrack <- GeneRegionTrack(txdb, 
                             name="Gene",
                             chromosome=chrs,
                             from = from ,
                             to = to )
  ht <- HighlightTrack(trackList=gtrack, 
                       start=as.numeric(region[5]), 
                       width=2, 
                       chromosome=chrs)
  pdf(file.path(fig.path,rbp,paste0(id,'.log.pdf')),width = 8, height = 8 )
  plotTracks(c(ht,grtrack,aligns),
               from = from ,
               to = to ,
               chromosome = chrs,
               sizes = c(0.4,0.6,rep(1,n)),
               #when error :Too many stacks to draw. ,
               #remenber expand grtrack size.
               extend.right = 50,
               extend.left = 50,
               type="coverage",
               showSampleNames=TRUE,
               main = id,
               cex.main = 1,
               transformation=function(x){log10(x+1)},
               labelPos="above")
    dev.off()
}

#===============================================================
sashimi.plot(id='ENSG00000100614:14:+:60277088:60293189:60299087',rbp ='CPSF6')
sashimi.plot(id='ENSG00000100614:14:+:60285088:60293189:60299087',rbp ='CPSF6')
sashimi.plot.log(id='ENSG00000100614:14:+:60285088:60293189:60299087',rbp ='CPSF6')
#===============================================================
##


python ./src/retrieve_seq.py /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./Analysis/PWMenrich/DAZAP1.ce.site ./Analysis/PWMenrich/de.fa ./Analysis/PWMenrich/de.txt