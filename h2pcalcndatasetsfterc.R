library(reshape2)

mappedERC <- function(ercmat, mapnames){
     map <- names(mapnames)
     names(map) <- mapnames
     multimap <- which(table(mapnames)>1)
     map[multimap] <- NA
     rownames(ercmat) <- map[rownames(ercmat)]
     colnames(ercmat) <- map[colnames(ercmat)]
     ercmat
}

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
t1 <- getERCasmat(merc, testset)
rm('merc')
mt1 <- melt(t1)
rm('t1')
mt1$value[!is.na(mt1$value)] <- 1
vdwyerc <- readRDS('../data/vertflywormyeast.ndatasetsfterc.rds')
mt2 <- melt(vdwyerc)
rm('vdwyerc')
mtc <- rbind(mt1,mt2)
rm('mt1')
rm('mt2')
tfinal <- acast(mtc, Var1 ~ Var2, sum, na.rm = T)
saveRDS(tfinal, file = '../data/all.ndatasetsfterc.rds')