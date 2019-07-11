args <- commandArgs(trailingOnly = T)

require(dplyr)
require(data.table)

source('../code/funcsERCfromProjections.R')

hsap.genes <- readLines('../data/hsap.genes')
genepresenceint <- readRDS('../data/genepresenceasinteger.rds')
intfterc <- readRDS('../data/allsum.fterc.rds')
pairmap <- readRDS('../data/pairmapmat.rds')
#merc.sym <- readRDS('../data/mamm63nt.trees.fterc.symm.rds')
zfg <- sort(readLines('../data/public/zincfinger.gene'))
org <- sort(readLines('../data/public/olfactoryreceptor.gene'))
orgzfg <- sort(c(org,zfg))
hsap.genes.useforsim <- setdiff(hsap.genes,orgzfg)

cplxf <- args[1]
cplxpairf <- args[2]
arrnum <- as.numeric(args[3])

cplx <- readRDS(cplxf)$cplx
rand.cplxsq <- fread(cplxpairf)[,c(1:2)]

nrandoms = 10
randcplxs <- generate.randcplxs(cplx, hsap.genes.useforsim, genepresenceint, nrandoms = nrandoms)
randcplxs.interc <- do.call('rbind',lapply(1:nrandoms, function(y){
     print(y)
     xv <- calc.betweencomplex.interc(rand.cplxsq,randcplxs[[y]])
     cbind(rand.cplxsq,xv) %>% mutate(dataset = ((arrnum-1)*nrandoms+y))
}))

outf <- gsub('.rds','',tail(strsplit(cplxf,'/',fixed = T)[[1]],1))
resdir <- paste0('../data/randoms/betweencomplex/',outf)
dir.create(resdir, showWarnings = F, recursive = T)
write.table(randcplxs.interc, file = paste0(resdir,'/runs',arrnum,'.tsv'),
            quote = F, sep = '\t', row.names = F)

