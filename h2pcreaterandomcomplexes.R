args <- commandArgs(trailingOnly = T)

require(dplyr)
require(data.table)

source('../code/funcsERCfromProjections.R')

hsap.genes <- readLines('../data/hsap.genes')
genepresenceint <- readRDS('../data/genepresenceasinteger.rds')
intfterc <- readRDS('../data/allsum.fterc.rds')
#pairmap <- readRDS('../data/pairmapmat.rds')
merc.sym <- readRDS('../data/mamm63nt.trees.fterc.symm.rds')
zfg <- sort(readLines('../data/public/zincfinger.gene'))
org <- sort(readLines('../data/public/olfactoryreceptor.gene'))
orgzfg <- sort(c(org,zfg))
hsap.genes.useforsim <- setdiff(hsap.genes,orgzfg)

cplxf <- args[1]
arrnum <- as.numeric(args[2])

cplx.u <- readRDS(cplxf)
nrandoms <- 500
randcplxs <- do.call('rbind',lapply(1:nrandoms, function(y){
     print(y)
     createRandomCplx.df(hsap.genes.useforsim, cplx.u$cplx, genepresenceint) %>% mutate(dataset = ((arrnum-1)*nrandoms+y))
}))
#randcplxs$pair <- with(randcplxs,calcpairid(V1,V2,pairmap))
# badpairids <- randcplxs$pair %in% badblastpairs$pair
# randcplxs$interc[badpairids] <- NA
# randcplxs$merc[badpairids] <- NA
randcplxs$interc <- with(randcplxs,calcpairid(V1,V2,intfterc))
randcplxs$merc <- with(randcplxs,calcpairid(V1,V2,merc.sym))

randcplxs.summ <- group_by(randcplxs, cplxid,dataset) %>%
     summarise(meanintscore = mean(interc, na.rm = T),
               meanmerc = mean(merc, na.rm = T))

#cplx <- '../data/public/corumcplxs.cleaned.rds'
outf <- gsub('.rds','',tail(strsplit(cplxf,'/',fixed = T)[[1]],1))
resdir <- paste0('../data/randoms/',outf)
dir.create(resdir, showWarnings = F, recursive = T)
write.table(randcplxs.summ, file = paste0(resdir,'/runs',arrnum,'.tsv'),
            quote = F, sep = '\t', row.names = F)
