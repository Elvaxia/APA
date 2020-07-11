#function for ALL
#=====================================================================
###sashimiplot
##=====================================================================
library(GenomicFeatures)
library(grid)
library(RColorBrewer)
library(Gviz)
options(ucscChromosomeNames=FALSE)
#load gtf
gtrack <- GenomeAxisTrack()
txdb <- loadDb('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/GRCh38.93.txdb')

#set colors for celltypes
ncolors <- length(unique(sample.info$Group))
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(ncolors)
names(my_colors) <- unique(sample.info$Group)

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
  align.list[s.name] <- s.align
}



path = '/mnt/raid61/Personal_data/xiaoxia/APA/Project/ALL_apa/Analysis/figures/sashimiplots/'
#plot for group
sashimi.plot <- function(id){
  region <- strsplit(id,':')[[1]]
  from = as.numeric(region[4])
  to = as.numeric(region[6])
  chrs <- region[2]
  aligns <- align.list[c(1,2,11,12)]
  grtrack <- GeneRegionTrack(txdb, 
                             name="Gene",
                             chromosome=chrs,
                             from = from ,
                             to = to )
  ht <- HighlightTrack(trackList=gtrack, 
                       start=as.numeric(region[5]), 
                       width=2, 
                       chromosome=chrs)
  pdf(paste0(path,id,'.pdf'),width = 8, height = 8 )
  plotTracks(c(ht,grtrack,aligns),
               from = from ,
               to = to ,
               chromosome = chrs,
               sizes = c(0.4,0.3,rep(1,4)),
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
#===============================================================


#===============================================================
##