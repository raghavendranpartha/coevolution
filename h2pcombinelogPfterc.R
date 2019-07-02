library(reshape2)

getERCasmat <- function(ercmat,genelist){
     #mat <- merc
     genv <- match(genelist,rownames(ercmat))
     print(genv)
     genord <- genelist[!is.na(genv)][order(genv[!is.na(genv)])]
     #out <- ercmat[genord,genord][lower.tri(ercmat[genord,genord], diag = F)]
     out <- ercmat[genord,genord]
     out[upper.tri(out)] <- 0
     out <- out+t(out)
     diag(out) <- 1
     out
}
merc <- readRDS('../data/mamm63nt.trees.fterc.rds')
testset <- readLines('../data/hsap.genes')
t1 <- -1*pnorm(getERCasmat(merc, testset), lower.tail = F, log.p = T)
rm('merc')
mt1 <- melt(t1)
rm('t1')
vdwyerc <- readRDS('../data/vertflywormyeast.logpfterc.rds')
mt2 <- melt(vdwyerc)
rm('vdwyerc')
mtc <- rbind(mt1,mt2)
rm('mt1')
rm('mt2')
tfinal <- acast(mtc, Var1 ~ Var2, sum, na.rm = T)
saveRDS(tfinal, file = '../data/allsum.logpfterc.rds')