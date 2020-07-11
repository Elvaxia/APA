#coregulation RBPs
library(STRINGdb)
library(data.table)
library(here)
library(dplyr)
load(here('Analysis/rawdata/rbp.DE.RData'))
rbp.DE.dt <- rbp.DE[[rbp]]

#load gene expression data
gene.path <- '/mnt/raid63/download/encode/KD-RNA-seq/hg19/expr_DESeq_combated/'
cpsf6.gene.exp <- fread('/mnt/raid63/download/encode/KD-RNA-seq/hg19/expr_DESeq_combated/ENCFF332PUX.tsv',sep='\t')
cpsf6.gene.exp <- 
cpsf6.de.gene <- cpsf6.gene.exp%>%filter(padj<0.1)%>%pull(id)
#map ensembl id back to gene name ,notice that gene expression used hg19
library(biomaRt)
grch37     <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                      path="/biomart/martservice",dataset="hsapiens_gene_ensembl")

ensemblids =gsub('[.][1-9]','',cpsf6.de.gene)
gene.names <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
      filters = 'ensembl_gene_id', 
      values = ensemblids, 
      mart = grch37)

Names <- gene.name[gene.name$ensembl_gene_id%in%ensemblids,]$external_gene_name

#PWMenriched RBPs
PWMenriched.RBPs <- fread(here('Analysis/PWMenrich/cpsf6.meme.res/centrimo.tsv'),header=T,fill=T)
colnames(PWMenriched.RBPs)[5] <- 'E_value'

selected.RBPs <- PWMenriched.RBPs%>% filter(E_value <1)%>%pull(motif_alt_id)%>%unique()
selected.RBPs <- c("CPSF6",selected.RBPs)
use.genes<- intersect(selected.RBPs,gene.names$external_gene_name)
use.genes <- c("CPSF6","HNRNPH2","RBM8A","LIN28A","SFPQ","SAMD4A")

string_db <- STRINGdb$new(version = "10", species = 9606, score_threshold = 400, input_directory = "")
use.genes <- as.data.frame(use.genes)
deg_mapped <- string_db$map(use.genes, "use.genes", removeUnmappedRows = TRUE )
hit <- deg_mapped$STRING_id
length(hit)
pdf(here('Analysis/figures/cpsf6.pro.interactive.pdf'))
string_db$plot_network(hit)
dev.off()

info <- string_db$get_interactions(deg_mapped$STRING_id)

# 采用igraph 进行聚类分析
clustersList <- string_db$get_clusters(deg_mapped$STRING_id)

string_db$plot_network(hit)

# 设置绘图参数
options(SweaveHooks=list(fig=function()
  par(mar=c(2.1, 0.1, 4.1, 2.1))))

# 绘制前4个聚类图
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}