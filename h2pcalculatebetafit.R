args <- commandArgs(trailingOnly=TRUE)
require(fitdistrplus)
ii <- as.numeric(args[1])
ercs <- c('/zfs1/nclark/rap119/erc/data/saves/mamm63nt.trees.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/vert39.trees.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/alldroso22.tre.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/worms17.trees.erc.rds',
          '/zfs1/nclark/rap119/erc/data/saves/all_Scer.tre.erc.rds')
minspthr <- c(10,10,5,5,5)
#for(ii in c(1:5)){
Sys.time()
erc <- readRDS(ercs[ii])
minsp <- minspthr[ii]
usepairs <- which(erc>=minsp, arr.ind = T)
nbrs <- 5
qqs = quantile(erc[usepairs], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
qqs[6] <- qqs[6]+1
betafit <- list()
for(jj in c(1:5)){
     print(jj)
     #ercdat <- erc[maxn>=qqs[jj] & maxn <qqs[jj+1]]
     ercdat <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
     ercdat <- ercdat[is.finite(ercdat)]
     ercdatscaled <- (ercdat+1.01)/2.02
     rm('ercdat')
     #rm('erc')
     betafit[[paste0(qqs[jj],'_',qqs[jj+1])]] <- fitdist(ercdatscaled, distr = 'beta',method = 'mle',lower = c(0, 0), keepdata = F)$est
     print(betafit)
}
#legend('topleft',legend = sapply(1:nbrs, function(x){paste(qqs[x:(x+1)], collapse = '_')}),col = cbPalette[1:5],lty = 1,lwd = 3)
fn <- gsub('.rds','',tail(strsplit(ercs[ii],'/')[[1]],1))
saveRDS(betafit, paste0('/zfs1/nclark/rap119/erc/data/saves/ercbetafit/',fn,'.betafits.rds'))
Sys.time()
#}