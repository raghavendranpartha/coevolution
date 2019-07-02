
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.fterc.rds')
werc <- readRDS('~/Documents/erc/data/saves/worms17.trees.fterc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')
mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')
testset <- readLines('../data/hsap.genes')
source('~/Documents/erc/code/funcsERCfromProjections.R')

t3 <- getERCasmat(mappedERC(derc, mapvecs[[1]]), testset)
t4 <- getERCasmat(mappedERC(werc, mapvecs[[2]]), testset)
t5 <- getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset)