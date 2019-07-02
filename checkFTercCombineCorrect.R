source('~/Documents/erc/code/funcsERCfromProjections.R')

genepres <- readRDS('../data/fulldatasetintegration_correct.rds')
genepresenceeach <- genepres$genepresenceeach
rm('genepres')
genepresenceeach[,1:5]
d1 <- which(colSums(genepresenceeach) == 1)
d2 <- which(colSums(genepresenceeach) == 2)
d3 <- which(colSums(genepresenceeach) == 3)
d4 <- which(colSums(genepresenceeach) == 4)
d5 <- which(colSums(genepresenceeach) == 5)
testset <- colnames(genepresenceeach)[c(d1[1:2],d2[1:2],d3[1:2],d4[2:3],d5[3:4])]

mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')

tall <- readRDS('../data/allsum.fterc.rds')

t1 <- getERCasmat(merc, testset)
t2 <- getERCasmat(verc, testset)
t3 <- getERCasmat(mappedERC(derc, mapvecs[[1]]), testset)
t4 <- getERCasmat(mappedERC(werc, mapvecs[[2]]), testset)
t5 <- getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset)

for(ii in seq(1,10,by=2)){
     ii <- 9
     print(genepresenceeach[,testset[ii:(ii+1)]])
     print(getERCasmat(tall, testset[ii:(ii+1)]))
     print(t1 <- getERCasmat(merc, testset[ii:(ii+1)]))
     print(t2 <- getERCasmat(verc, testset[ii:(ii+1)]))
     print(t3 <- getERCasmat(mappedERC(derc, mapvecs[[1]]), testset[ii:(ii+1)]))
     print(t4 <- getERCasmat(mappedERC(werc, mapvecs[[2]]), testset[ii:(ii+1)]))
     print(t5 <- getERCasmat(mappedERC(yerc, mapvecs[[3]]), testset[ii:(ii+1)]))
     
}


png('../figures/fterc_quantiles.png', width = 1000,height = 1000,res=120)
npairs <- 1000000
layout(matrix(c(5,0,4,0,3,0,2,0,1),3,3,byrow = T), widths = c(1,1,1),heights = c(1,1,1))
qqplot(rnorm(npairs),sample(yvec,size = npairs),ylim = c(-8,11),xlab = 'Normal theoretical quantiles', ylab = 'Sample quantiles', main = 'Yeast')
abline(a=0,b=1,lty=2,lwd=2,col='red')
qqplot(rnorm(npairs),sample(wercmap[lower.tri(wercmap)],size = npairs),ylim = c(-8,11),xlab = 'Normal theoretical quantiles', ylab = 'Sample quantiles', main = 'Worm')
abline(a=0,b=1,lty=2,lwd=2,col='red')
qqplot(rnorm(npairs),sample(merc[lower.tri(merc)],size = npairs),ylim = c(-8,11),xlab = 'Normal theoretical quantiles', ylab = 'Sample quantiles', main = 'Mammal')
abline(a=0,b=1,lty=2,lwd=2,col='red')
qqplot(rnorm(npairs),sample(dercmap[lower.tri(dercmap)],size = npairs),ylim = c(-8,11),xlab = 'Normal theoretical quantiles', ylab = 'Sample quantiles', main = 'Fly')
abline(a=0,b=1,lty=2,lwd=2,col='red')
qqplot(rnorm(npairs),sample(verc[lower.tri(verc)],size = npairs),ylim = c(-8,11),xlab = 'Normal theoretical quantiles', ylab = 'Sample quantiles', main = 'Vertebrate')
abline(a=0,b=1,lty=2,lwd=2,col='red')
dev.off()

