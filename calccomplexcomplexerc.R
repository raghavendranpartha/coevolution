require(dplyr)
require(data.table)

marc.cplx <- readRDS('../data/public/marcottecplxs.cleaned.rds')
length(marc.cplx$cplx)
min(sapply(marc.cplx$cplx, length))
max(sapply(marc.cplx$cplx, length))

marc.cplxsq <- as.data.frame(t(combn(1:length(marc.cplx$cplx),2)))
marc.cplxsq <- cbind(marc.cplxsq,calc.betweencomplex.interc.numpairs(marc.cplxsq,marc.cplx$cplx))
setnames(marc.cplxsq,c('complexind1','complexind2','meaninterc','numgenepairs'))
marc.cplxsq <- filter(marc.cplxsq, numgenepairs > 2, !is.na(meaninterc))

write.table(marc.cplxsq, file = '../data/public/betweencomplex/marcotte.betweencomplex.tsv', quote = F, sep = '\t', row.names = F)

corum.cplx <- readRDS('../data/public/corumcplxs.cleaned.rds')
corum.cplxsq <- as.data.frame(t(combn(1:length(corum.cplx$cplx),2)))
corum.cplxsq <- cbind(corum.cplxsq,calc.betweencomplex.interc.numpairs(corum.cplxsq,corum.cplx$cplx))
setnames(corum.cplxsq,c('complexind1','complexind2','meaninterc','numgenepairs'))
corum.cplxsq <- filter(corum.cplxsq, numgenepairs > 2, !is.na(meaninterc))

write.table(corum.cplxsq, file = '../data/public/betweencomplex/corum.betweencomplex.tsv', quote = F, sep = '\t', row.names = F)










genes <- hsap.genes.useforsim
genepresenceint <- readRDS('../data/genepresenceasinteger.rds')

marc.randcplxs <- generate.randcplxs(marc.cplx$cplx, hsap.genes.useforsim, genepresenceint, nrandoms = 100)
marc.rand.cplxsq <- marc.cplxsq[,c(1:2)]
marc.randcplxs.interc <- do.call('rbind',lapply(1:2, function(y){
     print(y)
     xv <- calc.betweencomplex.interc(marc.rand.cplxsq,marc.randcplxs[[y]])
     cbind(marc.rand.cplxsq,xv) %>% mutate(dataset = y)
}))


