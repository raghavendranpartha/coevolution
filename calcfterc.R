ercs <- c('~/Documents/erc/data/saves/mamm63nt.trees.erc.rds', '~/Documents/erc/data/saves/vert39.trees.erc.rds', '~/Documents/erc/data/saves/alldroso22.tre.erc.rds', '~/Documents/erc/data/saves/worms17.trees.erc.rds','~/Documents/erc/data/saves/all_Scer.tre.erc.rds')
ftercs <- c('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds', '~/Documents/erc/data/saves/vert39.trees.fterc.rds', '~/Documents/erc/data/saves/alldroso22.tre.fterc.rds', '~/Documents/erc/data/saves/worms17.trees.fterc.rds','~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')
titls <- c('mamm','vert','fly','worm','yeast')
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
par(mfrow=c(1,1))
for(jj in c(1:5)){
     Sys.time()
     #ii <- 1
     #tre <- readRDS(trs[ii])
     fterc <- readRDS(ftercs[jj])
     if(jj == 1){
          plot(density(fterc[lower.tri(fterc)], na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,0.5),xlim=c(-4,4),
               xlab = 'erc ft', ylab = 'density',main = 'FT ERC')
          #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
     }else{
          lines(density(fterc[lower.tri(fterc)], na.rm=T), col = cbPalette[jj], lwd = 2)
          #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
     }
     legend('topleft',legend = titls,
            col = cbPalette[1:5],lty = 1,lwd = 3,bty = 'n')
}
      
# ercft2 <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds')
# identical(ercft, ercft2)

     
titls <- c('mamm','vert','fly','worm','yeast')
cbPalette <-  c("#000000", "#E69F00","#0072B2", "#009E73","#D55E00", "#CC79A7","#F0E442","#56B4E9")
revinds <- function(x){
     out <- x
     out[,1] <- x[,2]
     out[,2] <- x[,1]
     out
}
#for(ii in c(1:5)){
#par(mfrow=c(2,2))
#for(ii in c(2:5)){
par(mfrow=c(1,1))
for(ii in c(1)){
     Sys.time()
     #ii <- 1
     #tre <- readRDS(trs[ii])
     erc <- readRDS(ercs[ii])
     minsp <- 10
     ercft <- erc
     goodinds <- which(erc >= minsp, arr.ind = T)
     ltinds <- sapply(1:nrow(goodinds), function(x){
          goodinds[x,1]>goodinds[x,2]
     })
     sum(ltinds)
     head(ercft[goodinds])
     revgoodinds <- revinds(goodinds)
     head(ercft[revgoodinds])
     ercft[lower.tri(ercft)] <- NA
     ercft[revgoodinds] <- atanh(erc[revgoodinds])*sqrt(ercft[goodinds]-3)
     head(ercft[revgoodinds])
     saveRDS(ercft, file = gsub('erc.rds','fterc.rds',ercs[ii]))
     #usepairs <- which(ercft>=minsp, arr.ind = T)
     nbrs <- 5
     qqs = quantile(ercft[goodinds], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     #rm('usepairs')
     qqs[6] <- qqs[6]+1
     
     for(jj in c(1:5)){
          #jj <- 3
          titl = titls[ii]
          #titl = ''
          print(jj)
          ercdat <- ercft[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          #nspdat <- t(ercft)[t(ercft)>=qqs[jj] & t(ercft) <qqs[jj+1]]
          ercdat <- ercdat[is.finite(ercdat)]
          #ercdatft <- atanh(ercdat)*sqrt(nspdat-3)
          #ercdatft <- ercdatft[is.finite(ercdatft)]
          #nspdat <- nspdat[is.finite(ercdatft)]
          if(jj == 1){
               plot(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2, ylim = c(0,0.5),xlim=c(-4,4),
                    xlab = 'erc ft', ylab = 'density',main = titl)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lty = 2, lwd = 2)
          }else{
               lines(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2)
               #lines(density(ercdatft, na.rm=T), col = cbPalette[jj], lwd = 2, lty = 2)
          }
          #rm('ercdat')
          rm('ercdat')
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),
            col = cbPalette[1:5],lty = 1,lwd = 3,bty = 'n')
     
}



#compare overall dist to sp binned
#lines(density(ercft[revgoodinds], na.rm = T), lwd = 2, col = 'magenta')
par(mfrow=c(2,2))
for(ii in c(2:5)){
     Sys.time()
     ii <- 1
     #tre <- readRDS(trs[ii])
     erc <- readRDS(ercs[ii])
     ercft <- readRDS(ftercs[ii])
     minsp <- 10
     goodinds <- which(erc >= minsp, arr.ind = T)
     nbrs <- 5
     qqs = quantile(ercft[goodinds], probs = seq(0,1,length.out = nbrs+1), na.rm = T)
     #rm('usepairs')
     qqs[6] <- qqs[6]+1
     plot(density(ercft[lower.tri(ercft)], na.rm=T),  lwd = 2, ylim = c(0,0.5),xlim=c(-4,4),col = 'red',lty=2,
          xlab = 'erc ft', ylab = 'density',main = titls[ii])
     for(jj in c(1:5)){
          print(jj)
          ercdat <- ercft[t(erc)>=qqs[jj] & t(erc) <qqs[jj+1]]
          
          ercdat <- ercdat[is.finite(ercdat)]
          print(length(ercdat))
          lines(density(ercdat, na.rm=T), col = cbPalette[jj], lwd = 2)
     }
     legend('topleft',legend = sapply(1:nbrs, function(x){paste0('[',paste(qqs[x:(x+1)], collapse = ','),')')}),
            col = cbPalette[1:5],lty = 1,lwd = 3,bty = 'n')
     
}

lines(density(ercft[lower.tri(ercft)], na.rm=T, adjust = 2),  lwd = 2, col = 'red')