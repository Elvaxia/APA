###sashimiplot
library(Gviz)
options(ucscChromosomeNames=FALSE)


align.list <- c()

ALL <- AlignmentsTrack('/mnt/data1/liyupeng/data/SecondHospital/BGI2/bam/CL100072036_L02_3.Aligned.sortedByCoord.out.bam',
                         isPaired = TRUE,
                         name="ALL",
                         fill="#5caf90")



Normal <- AlignmentsTrack('/mnt/data1/liyupeng/data/SecondHospital/BGI2/bam/CL100072036_L01_102.Aligned.sortedByCoord.out.bam',
                      isPaired = TRUE,
                      name="CTRL",
                      fill="#f7dc56")
gtrack <- GenomeAxisTrack()
library(GenomicFeatures)
txdb <- loadDb('/mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/GRCh38.93.txdb')
txTr <- GeneRegionTrack(txdb,name='Ensembl')
gtrack <- GenomeAxisTrack()
sashimi.plot <- function(id){
  #prepare bam files
  region <- strsplit(id,':')[[1]]
  from = as.numeric(region[4])
  to = as.numeric(region[6])
  chrs <- region[2]
  aligns <- list(Normal,ALL)
  n <- length(aligns)
  grtrack <- GeneRegionTrack(txdb,
                             name="Gene",
                             chromosome=chrs,
                             from = from ,
                             to = to )
  ht <- HighlightTrack(trackList=gtrack,
                       start=as.numeric(region[5]),
                       width=2,
                       chromosome=chrs)
  pdf(paste0('/mnt/raid61/Personal_data/xiaoxia/APA/Project/ALL_apa/Analysis/figures/sashimiplots/',id,'.pdf'),width = 8, height = 8 )
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
sashimi.plot('ENSG00000175643:16:+:11350654:11351387:11351762')
sashimi.plot('ENSG00000136536:2:+:159767343:159768515:159771027')


