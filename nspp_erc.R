library(RERconverge)
require(dplyr)
require(data.table)
require(fitdistrplus)

#plot density of erc values acc to # common branches 
trs <- c('~/Documents/erc/data/saves/mamm63nt.trees.rds', '~/Documents/erc/data/saves/vert39.trees.rds', '~/Documents/erc/data/saves/alldroso22.tre.rds' ,'~/Documents/erc/data/saves/worms17.trees.rds', '~/Documents/erc/data/saves/all_Scer.tre.rds')
ercs <- c('~/Documents/erc/data/saves/mamm63nt.trees.erc.rds', '~/Documents/erc/data/saves/vert39.trees.erc.rds', '~/Documents/erc/data/saves/alldroso22.tre.erc.rds', '~/Documents/erc/data/saves/worms17.trees.erc.rds','~/Documents/erc/data/saves/all_Scer.tre.erc.rds')
minspthr <- c(10,10,5,5,5)
titls <- c('mamm','vert','fly','worm','yeast')
ymaxs <- c(3.5,2.25,2.0,1.75,1.75)
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
#for(ii in c(1:5)){
par(mfrow=c(2,2))
for(ii in c(2:5)){
     Sys.time()
     ii <- 1
     #tre <- readRDS(trs[ii])
     erc <- readRDS(ercs[ii])
     minsp <- minspthr[ii]
     usepairs <- which(erc>=minsp, arr.ind = T)
     nbrs <- 5
     qqs = quantile(erc[usepairs], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     rm('usepairs')
     qqs[6] <- qqs[6]+1
     betafit <- list()
     for(jj in c(1:5)){
          #jj <- 3
          titl = titls[ii]
          print(jj)
          ercdat <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          ercdat <- ercdat[is.finite(ercdat)]
          if(jj == 1){
               plot(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,ymaxs[ii]),
                    xlab = 'erc', ylab = 'density',main = titl)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
          }else{
               lines(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
          }
          rm('ercdat')
          #rm('ercdatft')
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),col = cbPalette[1:5],lty = 1,lwd = 3, bty = 'n')
}

#density of fisher-transformed ercs
ymaxs <- c(3.5,2.25,2,1.75,1.75)
par(mfrow=c(2,2))
for(ii in c(2:5)){
     Sys.time()
     ii <- 1
     #tre <- readRDS(trs[ii])
     erc <- readRDS(ercs[ii])
     minsp <- minspthr[ii]
     usepairs <- which(erc>=minsp, arr.ind = T)
     nbrs <- 5
     qqs = quantile(erc[usepairs], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     rm('usepairs')
     qqs[6] <- qqs[6]+1
     betafit <- list()
     for(jj in c(1:5)){
          #jj <- 3
          titl = titls[ii]
          print(jj)
          ercdat <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          ercdat <- ercdat[is.finite(ercdat)]
          ercdatft <- atanh(ercdat)
          ercdatft <- ercdatft[is.finite(ercdatft)]
          if(jj == 1){
               plot(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,ymaxs[ii]),xlim=c(-2,2),
                    xlab = 'erc', ylab = 'density',main = titl)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
          }else{
               lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
          }
          #rm('ercdat')
          rm('ercdatft')
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),col = cbPalette[1:5],lty = 1,lwd = 3,bty = 'n')
}

# fit beta distribution to ercs;
for(ii in c(2:5)){
     Sys.time()
     ii <- 1
     #tre <- readRDS(trs[ii])
     erc <- readRDS(ercs[ii])
     #erc[1:6,1:6]
     minsp <- minspthr[ii]
     #nspp <- erc[upper.tri(erc)]
     usepairs <- which(erc>=minsp, arr.ind = T)
     nbrs <- 5
     qqs = quantile(erc[usepairs], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     rm('usepairs')
     qqs[6] <- qqs[6]+1
     betafit <- list()
     for(jj in c(1:5)){
          #jj <- 1
          #titl = titls[ii]
          print(jj)
          #ercdat <- erc[maxn>=qqs[jj] & maxn <qqs[jj+1]]
          ercdat <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          ercdat <- ercdat[is.finite(ercdat)]
          ercdatscaled <- (ercdat+1.01)/2.02
          #           if(jj == 1){
          #                plot(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,ymaxs[ii]),
          #                     xlab = 'erc', ylab = 'density',main = titl)
          #           }else{
          #                lines(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2)
          #           }
          rm('ercdat')
          #rm('erc')
          betafit[[paste0(qqs[jj],'_',qqs[jj+1])]] <- fitdist(ercdatscaled, distr = 'beta',method = 'mle',keepdata = F)$est
          print(betafit)
     }
     #legend('topleft',legend = sapply(1:nbrs, function(x){paste(qqs[x:(x+1)], collapse = '_')}),col = cbPalette[1:5],lty = 1,lwd = 3)
     fn <- gsub('.rds','',tail(strsplit(ercs[ii],'/')[[1]],1))
     saveRDS(betafit, paste0('~/Documents/erc/data/saves/ercbetafit/',fn,'.betafits.rds'))
     Sys.time()
}

cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
ymaxs <- c(6.5,4.5,4.0,2.75,3.25)
titls <- c('mamm','vert','fly','worm','yeast')
x <- seq(0,1,length.out = 201)
xsc <- 2.02*x-1.01
#ercbetafits <- grep('rds',dir('../data/ercbetafit/', include.dirs = FALSE),value = T)
ercbetafitsfs <- c('~/Documents/erc/data/saves/ercbetafit/mamm63nt.trees.erc.betafits.rds',
                   '~/Documents/erc/data/saves/ercbetafit/vert39.trees.erc.betafits.rds',
                   '~/Documents/erc/data/saves/ercbetafit/alldroso22.tre.erc.betafits.rds',
                   '~/Documents/erc/data/saves/ercbetafit/worms17.trees.erc.betafits.rds',
                   '~/Documents/erc/data/saves/ercbetafit/all_Scer.tre.erc.betafits.rds')
par(mfrow=c(2,2))
for(ii in 2:5){
     ii<-1
     bets <- readRDS(ercbetafitsfs[ii])
     for(jj in 1:5){
          if(jj == 1){
               plot(xsc,dbeta(x,shape1 = bets[[jj]][1], shape2 = bets[[jj]][2]), type = 'l', lwd = 3, col = cbPalette[jj], ylim = c(0,ymaxs[ii]),
                    xlab = 'erc',ylab = 'beta density', main = titls[ii])
          }else{
               lines(xsc,dbeta(x,shape1 = bets[[jj]][1], shape2 = bets[[jj]][2]), lty = 1, lwd = 3, col = cbPalette[jj])
          }
     }
     legend('topleft',paste0('[',gsub('_',',',names(bets)),')'),col=cbPalette[1:5],lty=1,lwd = 3,bty='n')
}






serc <- erc[1:6,1:6]

ybets <- readRDS('../data/ercbetafit/all_Scer.erc.betafits.rds')
xx <- list()
xx[[1]] <- c(2,3)
xx[[2]] <- c(4,2)


#load projdefault indiv

mmaxn <- mtre$report%*%t(mtre$report)
mmaxn[upper.tri(mmaxn, diag = T)] = NA
mmaxn[mmaxn < 10]=NA
vmaxn <- vtre$report%*%t(vtre$report)
vmaxn[upper.tri(vmaxn, diag = T)] = NA
vmaxn[vmaxn < 10]=NA
dmaxn <- dtre$report%*%t(dtre$report)
dmaxn[upper.tri(dmaxn, diag = T)] = NA
dmaxn[dmaxn < 5]=NA
wmaxn <- wtre$report%*%t(wtre$report)
wmaxn[upper.tri(wmaxn, diag = T)] = NA
wmaxn[wmaxn < 5]=NA
ymaxn <- ytre$report%*%t(ytre$report)
ymaxn[upper.tri(ymaxn, diag = T)] = NA
ymaxn[ymaxn < 5]=NA

cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
plotnsperc <- function(merc, mmaxn, titl = '',ymax = 3){
     #titl = 'Mamm'
     #ymax = 3
     merc[is.infinite(merc)] <- NA
     qqs = quantile(mmaxn, probs = seq(0,1,length.out = 6), na.rm = T)
     print(qqs)
     for(ii in c(1:5)){
          print(ii)
          if(ii == 1){
               #ii = 1
               #hist(merc[mmaxn>=qqs[ii]&mmaxn<=qqs[ii+1]], probability = T, nbreaks = 25, col = cbPalette[ii], na.rm = T)
               plot(density(merc[mmaxn>=qqs[ii]&mmaxn<=qqs[ii+1]], na.rm=T), col = cbPalette[ii], lwd = 2, ylim = c(0,ymax),
                    xlab = 'erc', ylab = 'density',main = titl)
          }else{
               #hist(merc[mmaxn>=qqs[ii]&mmaxn<=qqs[ii+1]], probability = T, nbreaks = 25, col = cbPalette[ii], na.rm = T, add = T)
               lines(density(merc[mmaxn>=qqs[ii]&mmaxn<=qqs[ii+1]], na.rm=T), col = cbPalette[ii], lwd = 3)
          }
          legend('topleft',legend = qqs[1:5],col = cbPalette[1:5],lty = 1,lwd = 3)
     }
}
plotnsperc(merc,mmaxn,'mamm',3.0)
plotnsperc(verc,vmaxn,'vert',2.2)
plotnsperc(derc,dmaxn,'dmel',2.0)
plotnsperc(werc,wmaxn,'worm',1.2)
plotnsperc(yerc,ymaxn,'scer',1.5)
