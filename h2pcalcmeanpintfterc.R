library(reshape2)

testset <- readLines('../data/hsap.genes')
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

mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')
testset <- readLines('../data/hsap.genes')

verc <- readRDS('../data/vert39.trees.fterc.rds')
t2 <- pnorm(getERCasmat(verc, testset), lower.tail = F, log.p = F)
rm('verc')
derc <- readRDS('../data/alldroso22.tre.fterc.rds')
t3 <- pnorm(getERCasmat(mappedERC(derc, mapvecs[[1]]), testset), lower.tail = F, log.p = F)
rm('derc')
werc <- readRDS('../data/worms17.trees.fterc.rds')
t4 <- pnorm(getERCasmat(mappedERC(werc, mapvecs[[2]]), testset), lower.tail = F, log.p = F)
rm('werc')
yerc <- readRDS('../data/all_Scer.tre.fterc.rds')
t5 <- pnorm(getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset), lower.tail = F, log.p = F)
rm('yerc')

tc <- rbind(melt(t2),melt(t3), melt(t4), melt(t5))
rm('t2')
rm('t3')
rm('t4')
rm('t5')

merc <- readRDS('../data/mamm63nt.trees.fterc.rds')
t1 <- pnorm(getERCasmat(merc, testset), lower.tail = F, log.p = F)
rm('merc')
mt1 <- melt(t1)
rm('t1')

mtc <- rbind(mt1,tc)
rm('mt1')
rm('tc')
tfinal <- acast(mtc, Var1 ~ Var2, mean, na.rm = T)
saveRDS(tfinal, file = '../data/allmeanp.fterc.rds')
