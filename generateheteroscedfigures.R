library(RERconverge)
require(dplyr)
require(data.table)
require(ggplot2)

source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
source('~/Documents/RERc/RERconverge/R/RcppExports.R')

mtre <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
vtre <- readRDS('~/Documents/erc/data/saves/vert39.trees.rds')
dtre <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')
wtre <- readRDS('~/Documents/erc/data/saves/worms17.trees.rds')
ytre <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.rds')

pdf('../figures/heterosced/mammal.pdf', width = 10,height = 6, family = 'FreeSans')
par(mfrow = c(1,2), mar = c(7,6,3,1), omi = c(0.6,0,0,0))
plotcomputeWeightsAllVar(mtre$paths, plot = T ,titl = 'Mammal')
dev.off()
pdf('../figures/heterosced/vertebrate.pdf', width = 10,height = 6, family = 'FreeSans')
par(mfrow = c(1,2), mar = c(7,6,3,1), omi = c(0.6,0,0,0))
plotcomputeWeightsAllVar(vtre$paths, plot = T ,titl = 'Vertebrate')
dev.off()
pdf('../figures/heterosced/Fly.pdf', width = 10,height = 6, family = 'FreeSans')
par(mfrow = c(1,2), mar = c(7,6,3,1), omi = c(0.6,0,0,0))
plotcomputeWeightsAllVar(dtre$paths, plot = T ,titl = 'Fly')
dev.off()
pdf('../figures/heterosced/worm.pdf', width = 10,height = 6, family = 'FreeSans')
par(mfrow = c(1,2), mar = c(7,6,3,1), omi = c(0.6,0,0,0))
plotcomputeWeightsAllVar(wtre$paths, plot = T ,titl = 'Worm')
dev.off()
pdf('../figures/heterosced/yeast.pdf', width = 10,height = 6, family = 'FreeSans')
par(mfrow = c(1,2), mar = c(7,6,3,1), omi = c(0.6,0,0,0))
plotcomputeWeightsAllVar(ytre$paths, plot = T ,titl = 'Yeast')
dev.off()
#plot lin unwt rers vs sqrt wt rers
plotcomputeWeightsAllVar=function (mat, plot = T,titl=''){
     set.seed(123)
     #mat <- ytre$paths
     transform = 'sqrt'
     
     plot=T
     nv=apply(mat, 2, mean,na.rm=T, trim=0.05)
     transform=match.arg(transform, choices = c("none", "sqrt", "log"))
     if (transform=="sqrt"){
          matsub=sqrt(mat)
          nvsub=sqrt(nv)
     }
     matr=naresidCPP(mat, model.matrix(~1+nv))
     matpred=fastLmPredictedMat(mat, model.matrix(~1+nv))
     mml=as.vector(mat)
     varl=as.vector(log(matr^2))
     matrsub=naresidCPP(matsub, model.matrix(~1+nvsub))
     matpredsub=fastLmPredictedMat(matsub, model.matrix(~1+nvsub))
     mmlsub=as.vector(matsub)
     varlsub=as.vector(log(matrsub^2))
     ii=which(!is.na(mml))
     mml=mml[ii]
     varl=varl[ii]
     iis=sample(length(mml), min(500000, length(mml)))
     mml=mml[iis]
     varl=varl[iis]
     mml2 = mml^2
     lmr2 <- lm(varl ~ mml + mml2)
     summary(lmr2)$r.sq
     l = lowess(mml,varl, f=0.7, iter = 2)
     f = approxfun(l, rule = 2)
     iisub=which(!is.na(mmlsub))
     mmlsub=mmlsub[ii]
     varlsub=varlsub[ii]
     iissub=sample(length(mmlsub), min(500000, length(mmlsub)))
     mmlsub=mmlsub[iissub]
     varlsub=varlsub[iissub]
     mml2sub = mmlsub^2
     lmr2sub <- lm(varlsub ~ mmlsub + mml2sub)
     summary(lmr2sub)$r.sq
     lsub = lowess(mmlsub,varlsub, f=0.7, iter = 2)
     fsub = approxfun(lsub, rule = 2)
     if (plot) {
          nbreaks=11
          qq=quantile(mml,seq(0,nbreaks,1)/nbreaks)
          qqdiff=diff(qq)
          breaks=qq[1:nbreaks]+qqdiff/2
          rr=quantile(mml, c(0.0001, 0.99))
          breaks=round(breaks,3)
          cutres<-cut(mml,breaks = breaks)
          cutres_tt=table(cutres)
          boxplot((varl)~ cutres, xlab = "", ylab = "log relative rate variance",
                  outline=F,  log="",xaxt='n',
                  cex.lab = 1.25,cex.axis = 1.25,lwd = 3)
          axis(1, at = 1:(nbreaks-1), labels = F)
          labels <- levels(cutres)
          ## Plot x axis labels at default tick marks
          text(1:(nbreaks-1), par("usr")[3] - 0.5, srt = 45, adj = 1,
               labels = labels, xpd = TRUE, cex = 1.25)
          title(paste0("A. ",titl," original relative rates"), adj = 0, font.main = 1, cex.main = 1.5)
          xx=(qq[1:nbreaks]+breaks)/2
          lines(1:length(xx), (f(qq[1:nbreaks])), lwd = 2, col = 2)
          #dev.off()
     }
     wr=1/exp(fsub(mmlsub))
     weights=(matrix(1/exp(fsub(matpredsub)), nrow = nrow(matsub)))
     if(plot){
          matrsub=naresidCPP(matsub, model.matrix(~1+nvsub), weights)
          varl2sub=(as.vector(log(matrsub^2))[iisub])[iissub]
          lmr2sub <- lm(varl2sub ~ mmlsub + mml2sub)
          summary(lmr2sub)$r.sq
          nbreaks=11
          qqsub=quantile(mmlsub,seq(0,nbreaks,1)/nbreaks)
          qqdiffsub=diff(qqsub)
          breakssub=qqsub[1:nbreaks]+qqdiffsub/2
          rrsub=quantile(mmlsub, c(0.0001, 0.99))
          breakssub=round(breakssub,3)
          cutressub<-cut(mmlsub,breaks = breakssub)
          l2sub = lowess(mmlsub,varl2sub, f=0.7, iter = 2)
          f2sub = approxfun(l2sub, rule = 2)
          boxplot((varl2sub)~ cutressub, xlab = "", ylab = "log relative rate variance",
                  outline=F,  log="",xaxt='n',
                  cex.lab = 1.25,cex.axis = 1.25,lwd = 3)
          axis(1, at = 1:(nbreaks-1), labels = F)
          labels <- levels(cutressub)
          text(1:(nbreaks-1), par("usr")[3] - 0.5, srt = 45, adj = 1,
               labels = labels, xpd = TRUE, cex = 1.25)
          lines(1:length(xx), (f2sub(qqsub[1:nbreaks])), lwd = 2, col = 2)
          title(paste0("B. ",titl," updated relative rates"), adj = 0, font.main = 1, cex.main = 1.5)
          #dev.off()
     }
     mtext('Average branch length',side = 1,outer = T,line = T, cex = 1.25)
}

