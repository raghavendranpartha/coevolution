##

#x <- matrix(sample(1:9,9), nrow = 3)

merc <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds')
verc <- readRDS('~/Documents/erc/data/saves/vert39.trees.fterc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.fterc.rds')
werc <- readRDS('~/Documents/erc/data/saves/worms17.trees.fterc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')
ercmatnames <- c('merc','verc','derc','werc','yerc')
titls <- c('Mammal','Vertebrate','Fly','Worm','Yeast')
minsp <- 10
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
#par(mfrow=c(2,2))
for(ii in c(1:5)){
     #ii <- 5
     pdf(paste0('../figures/ftercvnspp/',ercmatnames[ii],'.pdf'), width = 4, height=4, family = 'FreeSans')
     Sys.time()
     erc <- get(ercmatnames[ii])
     nspps <- erc[upper.tri(erc)]
     nbrs <- 5
     qqs = quantile(nspps[nspps >= 10], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     rm('nspps')
     qqs[6] <- qqs[6]+1
     for(jj in c(1:5)){
          titl = titls[ii]
          print(jj)
          ercdatft <-erc[lower.tri(erc) & t(erc) >= qqs[jj] & t(erc) < qqs[jj+1]]
          #erc[upper.tri(erc, diag = T)] <- NA
          #ercdatft <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          ercdatft <- ercdatft[is.finite(ercdatft)]
          if(jj == 1){
               plot(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,0.5),xlim=c(-3,3),
                    xlab = 'ERC', ylab = 'density',main='')
               title(titl, adj = 0, font.main = 1)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
          }else{
               lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
          }
          #rm('ercdat')
          rm('ercdatft')
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),
            col = cbPalette[1:5],pch = 19,bty = 'n')
     dev.off()
}

merc <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.erc.rds')
verc <- readRDS('~/Documents/erc/data/saves/vert39.trees.erc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.erc.rds')
werc <- readRDS('~/Documents/erc/data/saves/worms17.trees.erc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.erc.rds')
ercmatnames <- c('merc','verc','derc','werc','yerc')
titls <- c('Mammal','Vertebrate','Fly','Worm','Yeast')
ymaxs <- c(3.5,2.25,2,1.75,1.75)
minsp <- 10
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
#par(mfrow=c(2,2))
for(ii in c(1:5)){
     #ii <- 5
     pdf(paste0('../figures/ftercvnspp/',ercmatnames[ii],'_untransformed.pdf'), width = 4, height=4, family = 'FreeSans')
     Sys.time()
     erc <- get(ercmatnames[ii])
     nspps <- erc[upper.tri(erc)]
     nbrs <- 5
     qqs = quantile(nspps[nspps >= 10], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     rm('nspps')
     qqs[6] <- qqs[6]+1
     for(jj in c(1:5)){
          titl = titls[ii]
          print(jj)
          ercdatft <-erc[lower.tri(erc) & t(erc) >= qqs[jj] & t(erc) < qqs[jj+1]]
          #erc[upper.tri(erc, diag = T)] <- NA
          #ercdatft <- erc[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          ercdatft <- ercdatft[is.finite(ercdatft)]
          if(jj == 1){
               plot(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,ymaxs[ii]),xlim=c(-1,1),
                    xlab = 'Un-transformed ERC', ylab = 'density',main='')
               title(titl, adj = 0, font.main = 1)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
          }else{
               lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
          }
          #rm('ercdat')
          rm('ercdatft')
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),
            col = cbPalette[1:5],pch = 19,bty = 'n')
     dev.off()
}
