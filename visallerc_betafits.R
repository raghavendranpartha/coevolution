require(dplyr)
require(data.table)
require(ggplot2)
require(fitdistrplus)

#merc <- readRDS('~/Documents/erc/data/saves/mamm63.erc.rds')
#merc <- readRDS('~/Documents/erc/data/saves/vert39.erc.rds')
#merc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.erc.rds')
#merc <- readRDS('~/Documents/erc/data/saves/all_Scer.erc.rds')
#merc <- readRDS('~/Documents/erc/data/saves/worms17.erc.rds')

mdat <- merc[lower.tri(merc)]
#rm('merc')
mdat <- mdat[is.finite(mdat)] 
#mdatscaled <- (mdat+1)/2
#rm('mdat')
#betafit <- fitdist(mdatscaled, distr = 'beta',method = 'mle',lower = c(0, 0), keepdata = F)
#saveRDS(betafit, file = '../data/ercbetafit/wormercbetafit.rds')
#betafit$data <- NULL
#saveRDS(betafit, file = '../data/ercbetafit/vertercbetafit.rds')
#saveRDS(betafit, file = '../data/ercbetafit/dmelercbetafit.rds')
#saveRDS(betafit, file = '../data/ercbetafit/scerercbetafit.rds')


nps <- 1000000
mbet <- readRDS('../data/ercbetafit/mammalercbetafit.rds')
xbs <- rbeta(nps,mbet$estimate[1],mbet$estimate[2])*2-1
hist(mdat, breaks = 50, probability = T,main='Beta fit to mammal ERC', xlab = 'mammal erc')
lines(density(xbs),lwd=2,xlim=c(-1,1))
abline(v=0, lty = 2, lwd = 2)
rm('mdat')

nps <- 5000000
mbet <- readRDS('../data/ercbetafit/mammalercbetafit.rds')
xbs <- rbeta(nps,mbet$estimate[1],mbet$estimate[2])*2-1
plot(density(xbs),lwd=2,xlim=c(-1,1),'Beta fits to ERC', xlab = 'erc')
abline(v=0, lty = 2, lwd = 2)
#plot(pts, pbeta(pts,mbet$estimate[1],mbet$estimate[2]))


lins <- c('vert','dmel','worm','scer')
cols <- c('red','blue','magenta','green')
names(cols) <- lins
for(ii in lins){
     betafit <- readRDS(paste0('../data/ercbetafit/',ii,'ercbetafit.rds'))
     xbs <- rbeta(nps,betafit$estimate[1],betafit$estimate[2])*2-1
     lines(density(xbs),col=cols[ii], lwd = 2)
}
legend('topleft',c('mamm',lins),col = c('black',cols),lty = 1,lwd=2)


mambeta <- readRDS('../data/ercbetafit/mammalercbetafit.rds')



source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
source('~/Documents/RERc/RERconverge/R/RcppExports.R')
source('~/Documents/RERc/RERconverge/R/projection_coevo.R')
library(RERconverge)
source('~/Documents/erc/code/funcsERCfromProjections.R')

verc <- readRDS('~/Documents/erc/data/saves/vert39.trees.erc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.erc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.erc.rds')

ytre <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.rds')
ywts <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.weights.rds')

#colSums(dtreall$report)
#table(rowSums(dtreall$report))

#testset <- sample(names(dtreall$trees),5)
mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')
testset <- c('ELOVL5','ABCC1','ATP6V0A1','ABCB10','FOLH1','POLE2')
testsetynames <- mapvecs[[3]][c('ELOVL5','ABCC1','ATP6V0A1','ABCB10','FOLH1','POLE2')]
print(cbind(testset,testsetynames))


correlateERCGeneListTrees(ytre,mapvecs[[3]][testset], weighted = T, weights = ywts,min.sp = 5)
getERCasmat(yerc,mapvecs[[3]][testset])










source("https://bioconductor.org/biocLite.R")
biocLite('GOSemSim')

require(topGO)
require(GOSemSim)