source('~/Documents/erc/code/funcsERCfromProjections.R')

ftercs <- c('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds', '~/Documents/erc/data/saves/vert39.trees.fterc.rds', '~/Documents/erc/data/saves/alldroso22.tre.fterc.rds', '~/Documents/erc/data/saves/worms17.trees.fterc.rds','~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')

titls <- c('mamm','vert','fly','worm','yeast')
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")

ii <- 1
for(ii in c(1:5)){
     erc <- readRDS(ftercs[ii])
     qqs <- quantile(erc[lower.tri(erc)], probs = seq(0,1,length.out = 500), na.rm = T)
     print(tail(qqs))
}

sum(erc[lower.tri(erc)]>=5, na.rm = T)

mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')

for(ii in c(3:5)){
     ii <- 3
     erc <- readRDS(ftercs[ii])
     erc[1:10,1:10]
     map <- names(mapvecs[[ii-2]])
     names(map) <- mapvecs[[ii-2]]
     rownames(erc) <- map[rownames(erc)]
     colnames(erc) <- map[colnames(erc)]
     erc[1:10,1:10]
}


testset <- readLines('../data/hsap.genes')
#t1 <- getERCasmat(merc, testset)
t2 <- getERCasmat(verc, testset)
t3 <- getERCasmat(mappedERC(derc, mapvecs[[1]]), testset)
t4 <- getERCasmat(mappedERC(werc, mapvecs[[2]]), testset)
t5 <- getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset)

library(reshape2)

print(Sys.time())
#tcomb <- acast(rbind(melt(t1), melt(t2), melt(t3), melt(t4), melt(t5)), Var1 ~ Var2, sum, na.rm = T)
tcomb2 <- acast(rbind(melt(t3), melt(t4), melt(t5)), Var1 ~ Var2, sum, na.rm = T)
tcomb <- acast(rbind(melt(t2), melt(tcomb2)), Var1 ~ Var2, sum, na.rm = T)
sum(tcomb[lower.tri(tcomb)] == 0)
saveRDS(tcomb, file = '../data/vertflywormyeast.fterc.rds')
#tcomb2 <- acast(rbind(melt(t1), melt(t2)), Var1 ~ Var2, sum, na.rm = T)
print(Sys.time())

allfterc <- readRDS('../data/allsum.fterc.rds')
dim(allfterc)
