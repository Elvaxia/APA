#function for hema_mouse bulk data
#=====================================================================
###sashimiplot
##=====================================================================
library(GenomicFeatures)
library(grid)
library(RColorBrewer)
library(Gviz)
options(ucscChromosomeNames=FALSE)
#load gtf
generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  
  return(res)
}
gtrack <- GenomeAxisTrack()
txdb <- loadDb('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Mm/GRCm38.93.sqlite')

#set colors for celltypes
ncolors <- length(unique(sample.info$CellType))
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(ncolors)
my_colors <- c('#e26d43','#bb5f3f',"#377EB8",'#f0af55',"#8DD3C7",
               "#FED9A6" ,"#FFFFCC", "#FDCDAC" ,"#984EA3","#C0DA80",
               "#DBF1B9",'#8aab52',"#D1C2D2","#FB9A99",'#bc372d',
               '#bc372d','#f2a43a')
names(my_colors) <- unique(sample.info$CellType)

#prepare bam files
align.list <- list()
for(s in sample.info$sample){
  s.name <- s
  s.bam <- sample.info[sample==s,]$Bam
  color <- sample.info[sample==s,]$CellType
  s.align <- AlignmentsTrack(s.bam,
                             isPaired = TRUE,
                             name=s.name,
                             fill.coverage = my_colors[color] )
  align.list[s.name] <- s.align
}

path = '/mnt/raid61/Personal_data/xiaoxia/APA/Project/Hema_mouse/bulk_tianjin/Analysis/figures/sashimiplots/'
#plot for lineage
sashimi.plot <- function(id){
  region <- strsplit(id,':')[[1]]
  from = as.numeric(region[4])
  to = as.numeric(region[6])
  chrs <- region[2]
  grtrack <- GeneRegionTrack(txdb, 
                             name="Gene",
                             chromosome=chrs,
                             from = from ,
                             to = to )
  ht <- HighlightTrack(trackList=gtrack, 
                       start=as.numeric(region[5]), 
                       width=2, 
                       chromosome=chrs)
  lineages <- unique(sample.info$Lineage)
  lapply(lineages, function(x){
    celltypes <- sample.info[Lineage==x,]$sample
    aligns <- align.list[celltypes]
    n <- length(aligns)
    pdf(paste0(path,x,id,'.pdf'),width = 8, height = n+2 )
    plotTracks(c(ht,grtrack,aligns),
               from = from ,
               to = to ,
               chromosome = chrs,
               sizes = c(0.4,0.3,rep(1,n)),
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
  })
}
#===============================================================


#===============================================================
##