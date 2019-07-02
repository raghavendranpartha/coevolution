args <- commandArgs(trailingOnly = T)

library(GOSemSim)#this script is for version 1.26.0
require(dplyr)
require(data.table)

gogaf <- fread('../data/public/gogaf/goa_human.gaf', header = T)
gogaf.bygene <- group_by(gogaf, V3) %>% summarise(terms = paste(V5, collapse = ","))
mapvec <- readRDS('../data/public/gogaf/hsap.unmapgogaf.genes.rds')

infile <- args[1]
arrnum <- as.numeric(args[2])
#infile <- '../data/toperc/top.sumlogpgt15.bestbitscorepclt10.tsv'
indf <- fread(infile, header = T) %>%
     dplyr::select(c('g1','g2','pair','interc')) %>%
     na.omit()
ids1 <- indf$g1 %in% gogaf.bygene$V3
ids2 <- indf$g2 %in% gogaf.bygene$V3

indf$g1[!ids1] <- mapvec[indf$g1[!ids1]]
indf$g2[!ids2] <- mapvec[indf$g2[!ids2]]
indf <- indf %>% na.omit()

ntodo <- 10
st <- (arrnum-1)*ntodo+1
en <- arrnum*ntodo
onttype = 'MF'
gsm <- sapply(st:en, function(x){
     mgoSim(filter(gogaf, V3 == indf$g1[x])$V5, filter(gogaf, V3 == indf$g2[x])$V5, ont = onttype)
})
outf <- tail(strsplit(infile,'/')[[1]],1)
names(gsm) <- indf[st:en,]$pair
dir.create(paste0('../data/geneSim/',onttype,'/',outf), showWarnings = FALSE, recursive = T)
saveRDS(gsm, file = paste0('../data/geneSim/',onttype,'/',outf,'/pairs.',arrnum,'.rds'))
