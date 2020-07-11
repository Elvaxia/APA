#PWMemrich
setwd('/mnt/raid61/Personal_data/xiaoxia/APA/Project/RBP')
suppressPackageStartupMessages({
    library(here)
    library(data.table)
    library(dplyr)
    library(here)
    library(PWMEnrich)
    library(Biostrings)
})

load(here('Analysis/rawdata/rbp.DE.RData'))
 apa.factors <- c('CPSF6','CPSF7','CSTF2T','CSTF2','CELF1')
 rbp <- apa.factors[1]


load(here('Analysis/PWMenrich/All_PPMs.RData'))
# db <- gsub('@.+','',names(CISBP_PPMs))
# test.rbp <- intersect(db,names(files))
# rbp <- test.rbp[1]




#make fasta file containing 100bp DNA sequence around polyA sites.
#(DE.apa detected vs CE.apa occur in sample pairs,but not significantly different)
#DE.apa
all.site <- fread(paste0(here('Analysis/rawdata/'),rbp,'.filtered.RUD.txt'), select = 22)$id
DE.site <- rbp.DE.dt$id
write.table(DE.site, file= paste0(here('Analysis/PWMenrich/'),rbp,'.de.site'),
            col.names=F,row.names=F,quote=F)
CE.site <- setdiff(all.site, DE.site)
write.table(CE.site, file= paste0(here('Analysis/PWMenrich/'),rbp,'.ce.site'),
            col.names=F,row.names=F,quote=F)


system('/home/xiaoxia/anaconda3/bin/python ./src/retrieve_seq.py /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./Analysis/PWMenrich/CPSF6.de.site ./Analysis/PWMenrich/cpsf6.de.fa ./Analysis/PWMenrich/cpsf6.de.txt')

system('/home/xiaoxia/anaconda3/bin/python ./src/retrieve_seq.py /mnt/raid61/Personal_data/xiaoxia/APA/Project/Reference/Hs/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./Analysis/PWMenrich/CPSF6.ce.site ./Analysis/PWMenrich/cpsf6.ce.fa ./Analysis/PWMenrich/cpsf6.ce.txt')

library(Biostrings)
de.seq <- readDNAStringSet(here('Analysis/PWMenrich/cpsf6.de.fa'))
de.seq <- de.seq[!grepl("N", de.seq)]
ce.seq <- readDNAStringSet(here('Analysis/PWMenrich/cpsf6.ce.fa'))
bg.seq <- ce.seq[!grepl("N", ce.seq)]

PPMs <- c(RBPDB_PPMs, CISBP_PPMs, Ray_PPMs, RBNS_PPMs)

ids <- mapply(function(x) x[2], strsplit(names(PPMs), "@"))
ids[is.na(ids)] <- names(PPMs)[is.na(ids)]

PWMs = toPWM(PPMs, ids = ids, targets = mapply(function(x) x[1], strsplit(names(PPMs), "@")))
registerCoresPWMEnrich(10)

# bg.denovo = makeBackground(motifs = PWMs, organism = "hg19", type = "logn", quick = FALSE)
bg.custom.pvalue = makeBackground(motifs = PWMs, bg.seq = bg.seq, type = "pval",  bg.source = "apa")
save(bg.custom.pvalue, file = here('Analysis/PWMenrich/p.val.RData'))

bg.custom.logn = makeBackground(motifs = PWMs, bg.seq = bg.seq, type = "logn",  bg.source = "apa")


save(bg.custom.logn, file = here('Analysis/PWMenrich/logn.RData'))
#meme
system("/home/xiaoxia/software/meme/bin/centrimo ./Analysis/PWMenrich/cpsf6.de.fa ./Analysis/PWMenrich/Ray2013_rbp_Homo_sapiens.dna_encoded.meme --neg ./Analysis/PWMenrich/cpsf6.ce.fa --oc ./Analysis/PWMenrich/cpsf6.meme.res")