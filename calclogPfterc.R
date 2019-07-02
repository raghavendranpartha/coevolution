require(reshape2)

verc <- readRDS('~/Documents/erc/data/saves/vert39.trees.fterc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.fterc.rds')
werc <- readRDS('~/Documents/erc/data/saves/worms17.trees.fterc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')
mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')
testset <- readLines('../data/hsap.genes')

#t1 <- getERCasmat(merc, testset)
t2 <- -1*pnorm(getERCasmat(verc, testset), lower.tail = F, log.p = T)
rm('verc')
t3 <- -1*pnorm(getERCasmat(mappedERC(derc, mapvecs[[1]]), testset), lower.tail = F, log.p = T)
rm('derc')
t4 <- -1*pnorm(getERCasmat(mappedERC(werc, mapvecs[[2]]), testset), lower.tail = F, log.p = T)
rm('werc')
t5 <- -1*pnorm(getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset), lower.tail = F, log.p = T)
rm('yerc')

tc <- rbind(melt(t2),melt(t3), melt(t4), melt(t5))
rm('t2')
rm('t3')
rm('t4')
rm('t5')
tc$value[!is.na(tc$value)] <- 1
tcomb2sam <- acast(tc, Var1 ~ Var2, sum, na.rm = T)
saveRDS(tcomb2sam, file = '../data/vertflywormyeast.ndatasetsfterc.rds')


library(reshape2)

print(Sys.time())
tcomb2 <- acast(rbind(melt(t3), melt(t4), melt(t5)), Var1 ~ Var2, sum, na.rm = T)

tcomb <- acast(rbind(melt(t2), melt(tcomb2)), Var1 ~ Var2, sum, na.rm = T)
saveRDS(tcomb, file = '../data/vertflywormyeast.logpfterc.rds')
# pnorm(yercmap[5:15,5:15], lower.tail = F, log.p = T)
# yercmap[5:15,5:15]
# -1*pnorm(yercmap[5:15,5:15], lower.tail = F, log.p = T)

#run final combo with mammal on h2p
alllogpfterc <- readRDS('~/Documents/coevolution/data/allsum.logpfterc.rds')
diag(alllogpfterc) <- 1
alllogpfterc[alllogpfterc==0] <- NA
saveRDS(alllogpfterc, file = '~/Documents/coevolution/data/allsum.logpfterc.rds')
allsumfterc <- readRDS('~/Documents/coevolution/data/allsum.fterc.rds')
diag(allsumfterc) <- 1
allsumfterc[allsumfterc==0] <- NA
saveRDS(allsumfterc, file = '~/Documents/coevolution/data/allsum.fterc.rds')

alllogpfterc[1:7,1:7]
allsumfterc[1:7,1:7]
sum(allsumfterc == 0)

#calc which datasets have a pair


