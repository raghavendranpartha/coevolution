require(stringr)
require(dplyr)
require(data.table)
require(PRROC)

library(extrafont)
#font_import()
loadfonts(device="pdf")       #Register fonts for Windows bitmap output
fonts()

calcpairid <- function(g1,g2,mat){
     rowinds <- match(g1, rownames(mat))
     colinds <- match(g2, colnames(mat))
     mat[rowinds + nrow(mat) * (colinds - 1)]
}

ndataset.matched.controlpairs <- function(ndsvec, nPerLevelmult = 2){
     ndstab <- table(ndsvec)
     do.call('rbind',lapply(names(ndstab), function(x){
          sample_n(filter(bgdbig,ndataset == as.numeric(x)), ndstab[x]*nPerLevelmult)
     }))
}


commoninds <-function(A,B){
     dtA = data.table(A)
     dtB = data.table(cbind(B,as.numeric(0)),key=paste0("g",seq(len=ncol(B))))
     setnames(dtB, c(head(names(dtB), -1L), "found"))
     inds=which(dtB[dtA, list(found=ifelse(is.na(found), 0, 1))]==1)
     return(inds)
}

cleanedcomplex <- function(cplx, namesc = NULL){
     cplx.hsap.use <- sapply(cplx, function(x){
          intersect(x,hsap.genes.useforsim)
     })
     if(!is.null(namesc)){
          if(length(namesc) == length(cplx.hsap.use)){
               names(cplx.hsap.use) <- namesc
          }
     }
     cplx.hsap.use <- cplx.hsap.use[sapply(cplx.hsap.use, length)>2]
     cplx.hsap.use.df <- do.call('rbind',lapply(1:length(cplx.hsap.use), function(x){
          as.data.frame(t(combn(cplx.hsap.use[[x]],2))) %>% mutate(cplxid=x)
     }))
     cplx.hsap.use.df <- cplx.hsap.use.df %>%
          mutate(interc = calcpairid(V1,V2,intfterc)) %>%
          mutate(pair = calcpairid(V1,V2,pairmap)) %>%
          mutate(merc = calcpairid(V1,V2,merc.sym))
     list(cplx = cplx.hsap.use,
          df=cplx.hsap.use.df)
}

plotbox <- function(xx,yy,quantiles=T,breaksi = NULL,nbreaks=10,titl = "",xlab="",ylab=""){
     par(mar = c(10,5,3,1))
     if(quantiles){
          qq=quantile(xx,seq(0,nbreaks,1)/nbreaks,na.rm = T)
          qqdiff=diff(qq)
          breaks=qq[1:nbreaks]+qqdiff/2
          breaks=round(breaks,6)
          print(breaks)
     }else{
          breaks <- seq(min(xx,na.rm=T),max(xx,na.rm = T),length.out = nbreaks)
     }
     if(!is.null(breaksi)){
          breaks <- breaksi
          nbreaks <- length(breaks)
     }
     cutres<-cut(xx,breaks = breaks)
     cutres_tt <- table(cutres)
     boxplot((yy)~ cutres, xlab = xlab, ylab = ylab, outline=F,  log="", las=2, srt = 45, xaxt='n')
     axis(1, at=1:(nbreaks-1), labels = F)
     text(1:(nbreaks-1), y =-1, srt = 45, adj = 1,
          labels = levels(cutres), xpd = TRUE, cex = 1.25)
     text(1:(nbreaks-1), y =1000, srt = 30, adj = 0.5,
          labels = cutres_tt, xpd = TRUE, cex = 1.0)
     title(titl)
}

plotbox2 <- function(controls,xx,yy,breaksi = NULL, leglabsi = '',titl = "",xlab="",ylab="",ylimi = c(-3,8)){
     par(mar = c(6,5,3,1))
     if(is.null(breaksi)){
          breaks <- quantile(xx, probs = c(0.95,0.98,0.99,1))
          leglabs <- c('Control','Top 5-2%','Top 2-1%','Top 1%')
     }else{
          breaks = breaksi
          leglabs <- leglabsi
     }
     xx[xx<=breaks[1]] <- NA
     controlsxx <- -1000
     breaks <- c(controlsxx-1, breaks)
     print(breaks)
     nbreaks <- length(breaks)
     cutres<-cut(c(rep(controlsxx,length(controls)),xx),breaks = breaks, right = T,include.lowest = F)
     cutres_tt <- table(cutres)
     print(cutres_tt)
     boxplot(c(controls,yy)~ cutres, xlab = xlab, ylab = ylab, outline=F,  log="", las=2, srt = 45, xaxt='n',
             ylim = ylimi, cex.lab = 1.25, lwd = 2)
     abline(h=median(controls, na.rm = T), lty = 2, lwd = 2)
     axis(1, at=1:(nbreaks-1), labels = F)
     text(1:(nbreaks-1), y =ylimi[1]-1, srt = 45, adj = 1,
          labels = c('NA',levels(cutres)[-1]), xpd = TRUE, cex = 1.25)
     text(1:(nbreaks-1), y = ylimi[2], adj = 0.5,
          labels = leglabsi, xpd = TRUE, cex = 1.25)
     title(titl,adj=0,font.main=1,cex.main=1.5)
}
createRandomCplx <- function(genes, cplx){
     #genes <- hsap.genes
     #cplx <- cplx.hsap.use
     cplx.nunits <- sapply(cplx,length)
     ngenes <- sum(cplx.nunits)
     randngenes <- sample(genes, ngenes)
     randcplx <- split(randngenes, 
                       sample(rep(1:length(cplx), times = cplx.nunits),
                              size = ngenes))
     randcplx
}
createRandomCplx.df.naive <- function(genes, cplx){
     #genes <- hsap.genes
     #cplx <- cplx.hsap.use
     cplx.nunits <- sapply(cplx,length)
     ngenes <- sum(cplx.nunits)
     randngenes <- sample(genes, ngenes)
     randcplx <- split(randngenes, 
                       sample(rep(1:length(cplx), times = cplx.nunits),
                              size = ngenes))
     do.call('rbind',lapply(1:length(randcplx),function(x){
          as.data.frame(t(combn(randcplx[[x]],2))) %>% mutate(cplxid=x)
     }))
}
createRandomCplx.df <- function(genes, cplx, genepresenceint){
     #genes <- hsap.genes.useforsim
     #cplx <- marc.u$cplx
     randcplx <- sapply(1:length(cplx), function(x){
          #x <- 1
          genestosim <- genes[!genes %in% cplx[[x]]]
          genestosimint <- genepresenceint[genestosim]
          cplxgenesint <- genepresenceint[cplx[[x]]]
          tabcplxgenesint <- table(cplxgenesint)
          genesimforcplx <- sapply(1:length(tabcplxgenesint), function(y){
               sample(genestosim[genestosimint == as.numeric(names(tabcplxgenesint)[y])],
                            tabcplxgenesint[y])
          })
          #table(genepresenceint[unlist(genesimforcplx)])
          unlist(genesimforcplx, use.names = F)
     })
     do.call('rbind',lapply(1:length(randcplx),function(x){
          as.data.frame(t(combn(randcplx[[x]],2))) %>% mutate(cplxid=x)
     }))
}

getsingleERCasdf <- function(mat,genelist){
     melt(getERCasmat(mat, genelist)) %>%
          mutate(Var1=as.character(Var1)) %>%
          mutate(Var2=as.character(Var2)) %>%
          filter(Var1!=Var2)%>% 
          mutate(pair = calcpairid(Var1,Var2,pairmap)) %>%
          dplyr::select(c('pair','value'))
}

getERCasdf <- function(genelist){
     dflist <- lapply(c('merc','verc','dercmap','wercmap','yercmap'), function(x){
          df <- melt(getERCasmat(get(x), genelist)) %>%
               mutate(dataset = x)
     })
     do.call('rbind',dflist[sapply(dflist,nrow)>1]) %>%
          mutate(Var1=as.character(Var1)) %>%
          mutate(Var2=as.character(Var2)) %>%
          rowwise() %>% filter(Var1!=Var2)%>% 
          mutate(pair = paste(sort(c(Var1,Var2)), collapse = '_')) %>%
          dplyr::select(c(5,3,4)) %>% unique()
}

mappedERC <- function(ercmat, mapnames){
     map <- names(mapnames)
     names(map) <- mapnames
     multimap <- which(table(mapnames)>1)
     map[multimap] <- NA
     rownames(ercmat) <- map[rownames(ercmat)]
     colnames(ercmat) <- map[colnames(ercmat)]
     ercmat
}

getERCasmat <- function(ercmat,genelist){
     #mat <- merc
     genv <- match(genelist,rownames(ercmat))
     #print(genv)
     genord <- genelist[!is.na(genv)][order(genv[!is.na(genv)])]
     #out <- ercmat[genord,genord][lower.tri(ercmat[genord,genord], diag = F)]
     out <- ercmat[genord,genord]
     out[upper.tri(out)] <- 0
     out <- out+t(out)
     diag(out) <- NA
     out
}
getERCasmatlowertri <- function(ercmat,genelist){
     #mat <- merc
     genv <- match(genelist,rownames(ercmat))
     print(genv)
     genord <- genelist[!is.na(genv)][order(genv[!is.na(genv)])]
     #out <- ercmat[genord,genord][lower.tri(ercmat[genord,genord], diag = F)]
     out <- ercmat[genord,genord]
     out
}
runpr.intvsmamm <- function(true.erc.df, bgd.erc.df,plot=F,titl = '',xlimi = c(0,1),ypts,ymax){
     #setnames(marc.df.summ,c('cplxid','interc','merc'))
     #setnames(randcplxdfsumm,c('cplxid','dataset','interc','merc'))
     #true.erc.df <- marc.df.summ
     #bgd.erc.df <- randcplxdfsumm
     print(class(true.erc.df))
     true.erc.df$lbl <- 1
     bgd.erc.df$lbl <- 0
     colids <- c('lbl','merc','interc')
     nb.df <- rbind(dplyr::select(true.erc.df,colids), dplyr::select(bgd.erc.df,colids))
     nb.df$lbl <- factor(nb.df$lbl)
     probsi <- c(0,0.2,0.4,0.6,0.8,0.95,0.99)
     probsi <- c(seq(0,0.95,0.05),0.99)
     qs <- quantile(true.erc.df$interc, na.rm = T, probs = probsi)
     qsm <- quantile(true.erc.df$merc, na.rm = T, probs = probsi)
     truec <- bgdc <- rep(NA, length(qs))
     truecm <- bgdcm <- rep(NA, length(qsm))
     for(ii in 1:length(qs)){
          #ii <- qs[1]
          truec[ii] <- sum(true.erc.df$interc >= qs[ii], na.rm = T)
          bgdc[ii] <- sum(bgd.erc.df$interc >= qs[ii], na.rm = T)
          truecm[ii] <- sum(true.erc.df$merc >= qsm[ii], na.rm = T)
          bgdcm[ii] <- sum(bgd.erc.df$merc >= qsm[ii], na.rm = T)
          #min(true.erc.df$interc, na.rm = T)
          #sum(bgd.erc.df$interc >= min(true.erc.df$interc, na.rm = T), na.rm = T)
     }
     plot((1-probsi),(truec/bgdc)/(truec[1]/bgdc[1]), log='y',type = 'l', col = 'red',
          xlab = 'Recall',ylab = 'Fold change in Precision',cex.lab = 1.25, lwd = 2, yaxt = 'n',ylim = c(1,ymax))
     abline(h=1,lty=2,lwd=2)
     axis(2,at = ypts,labels = ypts)
     lines((1-probsi),(truecm/bgdcm)/(truecm[1]/bgdcm[1]), log='y',type = 'l', lwd = 2)
     legend('right', c('Integrated','Mammal'), lty = c(1), col = c('red','black'),
            bty = 'n',cex = 1.25, title = 'ERC', lwd = 2)
     title(titl, adj = 0, font.main = 1, cex.main = 1.5)
}
runroc.intvsmamm <- function(true.erc.df, bgd.erc.df,plot=F,titl = ''){
     true.erc.df$lbl <- 1
     bgd.erc.df$lbl <- 0
     colids <- c('lbl','merc','interc')
     nb.df <- rbind(dplyr::select(true.erc.df,colids), dplyr::select(bgd.erc.df,colids))
     nb.df$lbl <- factor(nb.df$lbl)
     fg <- filter(nb.df, lbl == 1)$interc
     bg <- filter(nb.df, lbl == 0)$interc
     fg_mammonly <- filter(nb.df, lbl == 1)$merc
     bg_mammonly <- filter(nb.df, lbl == 0)$merc
     pr <- roc.curve(scores.class0 = fg[!is.na(fg)], scores.class1 = bg[!is.na(bg)], curve = T)
     pr2 <- roc.curve(scores.class0 = fg_mammonly[!is.na(fg_mammonly)], scores.class1 = bg_mammonly[!is.na(bg_mammonly)], curve = T)
     if(plot){
          plot(pr$curve[,1],pr$curve[,2],xlim = c(0,1), type = 'l',col='red',xlab='FPR',ylab='TPR',lwd = 2,main = '',
               cex.lab = 1.25)
          title(titl, font.main = 1, cex.main = 1.5, adj = 0)
          lines(pr2$curve[,1],pr2$curve[,2],xlim = c(0,1), lwd = 2)
          abline(a = 0,b=1,lty = 2)
          legend('bottomright', c(paste0('Integrated: AUC = ',round(pr$auc, 2)),paste0('Mammal: AUC = ',round(pr2$auc, 2))),
                 title = 'ERC', pch = 19,col = c('red','black'),bty='n', title.adj = 0.5,cex=1.25)
     }
     return(list(int.auc = pr$auc,mamm.auc = pr2$auc))
}
runroc.string.intvsmamm <- function(true.erc.df, bgd.erc.df, levels = NULL,plot=F,titl = '',
                                    nPointsforcurve = NULL){
     #titl = 'STRING: Any'
     #true.erc.df <- string.type.use.forall
     #bgd.erc.df <- sample_n(onembgdpairs, size = 500000)
     #levels = c(930,951,960)
     true.erc.df$lbl <- 1
     bgd.erc.df$lbl <- 0
     bgd.erc.df$score <- -1000
     colids <- c('lbl','merc','interc','score')
     nb.df <- rbind(dplyr::select(true.erc.df,colids), dplyr::select(bgd.erc.df,colids))
     nb.df$lbl <- factor(nb.df$lbl)
     for(ii in 1:length(levels)){
          assign(paste0('fg',ii), filter(nb.df, lbl == 1, score > levels[ii])$interc)
          assign(paste0('bg',ii), filter(nb.df, lbl == 1)$interc)
          assign(paste0('fg_mammonly',ii), filter(nb.df, lbl == 1, score > levels[ii])$merc)
          assign(paste0('bg_mammonly',ii), filter(nb.df, lbl == 1)$merc)
          assign(paste0('rocint',ii),
                 roc.curve(scores.class0 = get(paste0('fg',ii))[!is.na(get(paste0('fg',ii)))], 
                           scores.class1 = get(paste0('bg',ii))[!is.na(get(paste0('bg',ii)))], curve = T))
          assign(paste0('rocmamm',ii),
                 roc.curve(scores.class0 = get(paste0('fg_mammonly',ii))[!is.na(get(paste0('fg_mammonly',ii)))], 
                           scores.class1 = get(paste0('bg_mammonly',ii))[!is.na(get(paste0('bg_mammonly',ii)))], curve = T))
     }
     ltys <- c(1,3,4)
     
     if(plot){
          cols <- c('blue','black','red')
          if(is.null(nPointsforcurve)){
               plot(rocint1$curve[,1],rocint1$curve[,2],xlim = c(0,1), type='l',lwd=2,col=cols[1],xlab='FPR',ylab='TPR',
                    cex.lab = 1.25, cex.axis = 1.25, main ='')
               title(titl, adj = 0, font.main = 1,cex.main=1.25)
               lines(rocmamm1$curve[,1],rocmamm1$curve[,2],xlim = c(0,1), lwd = 2, col = cols[1], lty = 5)
               for(ii in 2:length(levels)){
                    lines(get(paste0('rocint',ii))$curve[,1],get(paste0('rocint',ii))$curve[,2],
                          col=cols[ii], lwd = 2)
                    lines(get(paste0('rocmamm',ii))$curve[,1],get(paste0('rocmamm',ii))$curve[,2],
                          col=cols[ii], lwd = 2, lty = 5)
               }
               abline(a = 0,b=1,lty = 3, lwd = 2)
               legend('right', c('Top 1%','Top 2%','Top 5%'), col = rev(cols), pch=19, bty = 'n', cex = 1.25, title = 'Dataset')
               legend('bottomright', c('Integrated','Mammal'), lty = c(1,5), bty = 'n',cex = 1.25, title = 'ERC', lwd = 2)
          }else{
               if(!is.null(nPointsforcurve)){
                    inds1 <- sort(sample(1:length(rocint1$curve[,1]), size = nPointsforcurve))
                    inds2 <- sort(sample(1:length(rocmamm1$curve[,1]), size = nPointsforcurve))
                    plot(rocint1$curve[inds1,1],rocint1$curve[inds1,2],xlim = c(0,1), type='l',lwd=2,col=cols[1],xlab='FPR',ylab='TPR',
                         cex.lab = 1.25, cex.axis = 1.25, main ='')
                    title(titl, adj = 0, font.main = 1,cex.main=1.25)
                    lines(rocmamm1$curve[inds2,1],rocmamm1$curve[inds2,2],xlim = c(0,1), lwd = 2, col = cols[1], lty = 5)
                    for(ii in 2:length(levels)){
                         inds1 <- sort(sample(1:length(get(paste0('rocint',ii))$curve[,1]), size = nPointsforcurve))
                         inds2 <- sort(sample(1:length(get(paste0('rocmamm',ii))$curve[,1]), size = nPointsforcurve))
                         lines(get(paste0('rocint',ii))$curve[inds1,1],get(paste0('rocint',ii))$curve[inds1,2],
                               col=cols[ii], lwd = 2)
                         lines(get(paste0('rocmamm',ii))$curve[inds2,1],get(paste0('rocmamm',ii))$curve[inds2,2],
                               col=cols[ii], lwd = 2, lty = 5)
                    }
                    abline(a = 0,b=1,lty = 3, lwd = 2)
                    legend('right', c('Top 1%','Top 2%','Top 5%'), col = rev(cols), pch=19, bty = 'n', cex = 1.25, title = 'Dataset')
                    legend('bottomright', c('Integrated','Mammal'), lty = c(1,5), bty = 'n',cex = 1.25, title = 'ERC', lwd = 2)
               }
          }
          
     }
     return(c(lapply(paste0('rocint',1:length(levels)), function(x) get(x)$auc),
              lapply(paste0('rocmamm',1:length(levels)), function(x) get(x)$auc)))
}

runprecrec.intvsmamm <- function(true.erc.df, bgd.erc.df,plot=F,titl = ''){
     nreps <- 10
     true.erc.df <- corumcplx.use.df
     bgd.erc.df <- randcplxs[!badrandpairs,]
     titl = 'corum'
     scorelabels1 <- lapply(1:nreps, function(x){
          sc <- c(true.erc.df$interc,
                  filter(bgd.erc.df, dataset %% nreps == x)$interc)
          lb <- c(true.erc.df$lbl,
                  filter(bgd.erc.df, dataset %% nreps == x)$lbl)
          list(scores = sc, labels = lb)
     })
     scorelabels2 <- lapply(1:nreps, function(x){
          sc <- c(true.erc.df$merc,
                  filter(bgd.erc.df, dataset %% nreps == x)$merc)
          lb <- c(true.erc.df$lbl,
                  filter(bgd.erc.df, dataset %% nreps == x)$lbl)
          list(scores = sc, labels = lb)
     })
     scores <- c(sapply(scorelabels1,'[[',1),sapply(scorelabels2,'[[',1))
     labels <- join_labels(sapply(scorelabels1,'[[',2),sapply(scorelabels2,'[[',2))
     mets <- rep(c('Integrated','Mammal only'), length = nreps)
     dsids <- rep(1:nreps, each = 2)
     nbcorummmdat <- mmdata(scores = c(scores1a,scores1b,scores2a,scores2b),
                            labels = c(labels1a,labels1b,labels2a,labels2b),
                            modnames = mets,
                            dsids = dsids)
     mmcurves <- evalmod(nbcorummmdat,cb_alpha = 0.05)
     autoplot(mmcurves, "PRC", show_cb = TRUE )+coord_cartesian(xlim = c(0,0.25))+
          labs(fill="Method")+
          theme_bw()+ggtitle('Naive Bayes Prediction: Corum Complexes')+
          theme(axis.text=element_text(size = 18, face='bold'),
                axis.title=element_text(size = 20, face='bold'),
                plot.title=element_text(size = 22, face='bold'),
                strip.text = element_text(size = 18, face = 'bold'))
     
     colids <- c('lbl','merc','interc')
     nb.df <- rbind(dplyr::select(true.erc.df,colids), dplyr::select(bgd.erc.df,colids))
     nb.df$lbl <- factor(nb.df$lbl)
}
     


runNaiveBayes <- function(true.erc.df, bgd.erc.df,plot=F){
     #true.erc.df <- marc.erc.df
     #bgd.erc.df <- rand.cplx.erc.df
     colids <- c('lbl','merc','verc','dercmap','wercmap','yercmap')
     nb.df <- rbind(dplyr::select(true.erc.df,colids), dplyr::select(bgd.erc.df,colids))
     nb.df$lbl <- factor(nb.df$lbl)
     Naive_Bayes_Model=naiveBayes(lbl ~., data=nb.df)
     Naive_Bayes_Model_mammonly=naiveBayes(lbl ~ merc, data=nb.df)
     NB_Predictions=predict(Naive_Bayes_Model,nb.df[,c(-1)],type='raw')
     NB_Predictions_mammonly=predict(Naive_Bayes_Model_mammonly,nb.df[,'merc'],type='raw')
     nb.df$nbp <- NB_Predictions[,2]
     nb.df$nbp_mammonly <- NB_Predictions_mammonly[,2]
     
     fg <- filter(nb.df, lbl == 1)$nbp
     bg <- filter(nb.df, lbl == 0)$nbp
     
     fg_mammonly <- filter(nb.df, lbl == 1)$nbp_mammonly
     bg_mammonly <- filter(nb.df, lbl == 0)$nbp_mammonly
     # ROC Curve    
#      roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
#      plot(roc)
#      
     # PR Curve
     pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
     pr2 <- pr.curve(scores.class0 = fg_mammonly, scores.class1 = bg_mammonly, curve = T)
     #pr <- pr.curve(scores.class0 = marc.nb.df$nbp, weights.class0 = as.numeric(marc.nb.df$lbl), curve = T)
     #plot(pr, xlim = c(0,0.1))
     #plot(pr2, xlim = c(0,0.1), add = T)
     if(plot){
          plot(pr$curve[,1],pr$curve[,2],xlim = c(0,1), type = 'l',col='red', lwd = 2, xlab = 'Recall',ylab='Precision',main = 'Naive Bayes Prediction')
          lines(pr2$curve[,1],pr2$curve[,2],xlim = c(0,1), lwd = 2)
          legend('topright', c('Integrated','Mammal only'),lwd = 2,col = c('red','black'),bty='n')
     }
     
     
     return(list(int.score = nb.df$nbp,
                 mamm.score = nb.df$nbp_mammonly,
                 lbl = nb.df$lbl))
     
}






correlateERCGeneListTrees=function(treesObj, genelist, plot = T, cutoff=NULL, transform="sqrt", weighted=T, weights = NULL, useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05, xlims = NULL){
     
#      totalpairs <- choose(treesObj$numTrees,2)
#      treesObj = dtreall
#      genelist <- testset
#      cutoff<-NULL
#      transform="sqrt"
#      weighted=T
#      weights = NULL
#      useSpecies=NULL
#      min.sp=5
#      scale=T
#      plot=T
#      doOnly=NULL
#      maxT=NULL
#      scaleForPproj=F
#      mean.trim=0.05
     
     
     totalpairs <- choose(length(genelist),2)
     
     if(is.null(cutoff)){
          cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
          message(paste("cutoff is set to", cutoff))
     }
     if (weighted){
          if(is.null(weights)){
               message("Weights can not be NULL if weighted = T; computing weights")
               weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
          }
          residfunc=fastLmResidMatWeighted
     }else{
          residfunc=fastLmResidMat
     }
     
     if (is.null(useSpecies)){
          useSpecies=treesObj$masterTree$tip.label
          #mappedEdges=trees$mappedEdges
     }
     if(is.null(maxT)){
          maxT=treesObj$numTrees
     }
     if(transform!="none"){
          transform=match.arg(transform,c("sqrt", "log"))
          transform=get(transform)
     }else{
          transform=NULL
     }
     #cm is the names of species that are included in useSpecies and the master tree
     cm=intersect(treesObj$masterTree$tip.label, useSpecies)
     sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
     if (length(sp.miss) > 0) {
          message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                       collapse = ",")))
          
     }
     
     corout=matrix(NA,treesObj$numTrees,treesObj$numTrees)
     diag(corout) = 1
     
     geneis <- sort(match(genelist, names(treesObj$trees)))
     #genejs <- match(genelist, names(treesObj$trees))
     
     maxn2=treesObj$report%*%t(treesObj$report)
     maxn2[upper.tri(maxn2, diag = T)] = NA
     corout[which(maxn2 == 0)] <- Inf
     
     
     #      if(is.null(doOnly)){
     #           doOnly=1
     #      }else{
     #           maxT=1
     #      }
     pairsdone <- 0
     print(maxT)
     for (i in geneis){
          for(j in geneis){
               if(i < j){
                    next
               }
               #print(paste0(i,'',j))
               #print(j)
               #i <- geneis[2]
               #j <- geneis[1]
               #i <- 160
               #j <- 67
               if(is.na(corout[i,j])){
                    #tree1=treesObj$trees[[i]]
                    #get the common species, prune and unroot
                    both=intersect(intersect(treesObj$trees[[i]]$tip, treesObj$trees[[j]]$tip), cm)
                    nb=length(both)
                    if(nb>2){
                         tree1=unroot(pruneTree(treesObj$master,both))
                    }
                    #find all the genes that contain all of the species in tree1
                    allreport=treesObj$report[,both]
                    if(nb>1){
                         ss=rowSums(allreport)
                    }else{
                         ss=allreport
                    }
                    iiboth=which(ss==length(both))
                    ai2=which(maxn2[iiboth, iiboth]==nb, arr.ind = T)
                    
                    if(length(both)<min.sp){
                         corout[matrix(iiboth[ai2],nrow=nrow(ai2))]=Inf
                         ai2rev <- ai2
                         ai2rev[,1] <- ai2[,2]
                         ai2rev[,2] <- ai2[,1]
                         corout[matrix(iiboth[ai2rev],nrow=nrow(ai2rev))]=1
                         pairsdone <- pairsdone+nrow(ai2)
                         #message(paste0(pairsdone," done out of ",totalpairs))
                         next
                    }    
                    
                    
                    if(T){
                         
                         
                         ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
                         
                         ii= treesObj$matIndex[ee[, c(2,1)]]
                         
                         allbranch=treesObj$paths[iiboth,ii]
                         if(weighted){
                              allbranchw=weights[iiboth,ii]
                         }
                         if(scaleForPproj){
                              nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
                         }else{
                              nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
                         }
                         
                         iibad=which(allbranch<cutoff)
                         #don't scale
                         #allbranch=scaleMat(allbranch)
                         if(!is.null(transform)){
                              nv=transform(nv)
                              allbranch=transform(allbranch)
                         }
                         allbranch[iibad]=NA
                         
                         if(!scale){
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                                   
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv), allbranchw)
                                   
                              }
                         }else{
                              
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv),allbranchw)
                              }
                              
                              proj=scale(proj, center = F)
                              
                         }
                         namesai2 <- as.data.frame(cbind(names(iiboth)[ai2[,1]],
                                           names(iiboth)[ai2[,2]]), stringsAsFactors = F)%>%
                              mutate(m=1:n())
                         namesai2ingenelist <- filter(namesai2,
                                                      V1 %in% genelist, V2 %in% genelist)
                         #we have the projection
                         
                         #for(m in 1:nrow(ai2)){
                         for(nn in 1:nrow(namesai2ingenelist)){
                              #nn = 1
                              #print(m)
                              m <- namesai2ingenelist$m[nn]
                              k=sort(ai2[m,])[1]
                              l=sort(ai2[m,])[2]
                              #tmpcor=cor(proj[k,], proj[l,], method = 'p', use = 'pair')
                              re1 <- proj[k,]
                              re2 <- proj[l,]
                              goodinds <- !is.na(re1)&!is.na(re2)
                              if(sum(!is.na(proj[k,])&!is.na(proj[l,])) >= min.sp){
                                   cres=cor.test(win(proj[k,],3), win(proj[l,],3), method='pearson', exact=F)
                                   corout[iiboth[l], iiboth[k]]=cres$est
                                   if(plot){
                                        crln <- corout[iiboth[l], iiboth[k]]
                                        crln.nowin <- cor.test(proj[k,], proj[l,], method='pearson', exact=F)$est
                                        pointnames <- colnames(allbranch)[goodinds]
                                        pointnames[is.na(pointnames)] <- ''
                                        par(mfrow = c(1,1))
                                        if(is.null(xlims)){
                                             xlims <- range(win(re1[goodinds],3))
                                        }
                                        plot(win(re1[goodinds],3),win(re2[goodinds],3), main = with(namesai2ingenelist,
                                                                                      paste0(V1[nn],'_',V2[nn],'\nerc = ',round(crln,3))),
                                             xlab = namesai2ingenelist$V2[nn], ylab = namesai2ingenelist$V1[nn], xlim = xlims)
                                        text(win(re1[goodinds],3),win(re2[goodinds],3), labels = pointnames,cex = 0.8)
#                                         plot(re1[goodinds],re2[goodinds], main = with(namesai2ingenelist,paste0(V1[nn],'_',V2[nn],'\nerc = ',round(crln.nowin,3))),
#                                              xlab = namesai2ingenelist$V2[nn], ylab = namesai2ingenelist$V1[nn])
#                                         text(re1[goodinds],re2[goodinds], labels = pointnames)
                                        #legend('topleft', legend = c(paste0('erc = ',round(crln,3))))
                                   }
                              }else{
                                   corout[iiboth[l], iiboth[k]]=Inf
                              }
                              corout[iiboth[k], iiboth[l]]=sum(!is.na(proj[k,])&!is.na(proj[l,]))
                         }
                         pairsdone <- pairsdone+nrow(ai2)
                         #message(paste0(pairsdone," done out of ",totalpairs))
                    }
                    
               }
               
          }
     }
     rownames(corout) <- names(treesObj$trees)
     colnames(corout) <- names(treesObj$trees)
     corout[is.infinite(corout)] <- NA
     corout[geneis,geneis]
}
correlateERCTreesinProjection <- function(treesObj, dfprojns = NULL, cutoff=NULL, transform="sqrt", weighted=T, weights = NULL, useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05){
     
#      treesObj = dtre
#      cutoff=NULL
#      transform="sqrt"
#      weighted=T
#      weights = weights.sqrt
#      useSpecies=NULL
#      min.sp=3
#      scale=T
#      doOnly=NULL
#      maxT=NULL
#      scaleForPproj=F
#      mean.trim=0.05
#      dfprojns = readRDS('../data/saves/onekdroso20.tre.projns.rds')
#      genepairsprojn = readRDS('../data/saves/onekdroso20.tre.projnGenepair.rds')
     
     if(is.null(dfprojns)){
          stop('Projection objects not supplied')
     }
     if(is.null(cutoff)){
          cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
          message(paste("cutoff is set to", cutoff))
     }
     #cutoff = 0
     if (weighted){
          #weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
          if(is.null(weights)){
               message("Weights can not be NULL if weighted = T; computing weights")
               weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
          }
          residfunc=fastLmResidMatWeighted
     }else{
          residfunc=fastLmResidMat
     }
     # residfunc=naresid
     
     if (is.null(useSpecies)){
          useSpecies=treesObj$masterTree$tip.label
          #mappedEdges=trees$mappedEdges
     }
     if(is.null(maxT)){
          maxT=treesObj$numTrees
     }
     if(transform!="none"){
          transform=match.arg(transform,c("sqrt", "log"))
          transform=get(transform)
     }else{
          transform=NULL
     }
     
     
     
     #cm is the names of species that are included in useSpecies and the master tree
     cm=intersect(treesObj$masterTree$tip.label, useSpecies)
     sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
     if (length(sp.miss) > 0) {
          message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                       collapse = ",")))
          
     }
     
     corout=matrix(0,nrow(treesObj$paths),nrow(treesObj$paths))
     diag(corout)=1
     #corout[lower.tri(corout, diag = F)] <- 0
     pairsdone <- 0
     totalpairs <- 'N/A'
     
     #maxn=rowSums(treesObj$report[,cm])
     
     maxn2=treesObj$report%*%t(treesObj$report)
     maxn2[upper.tri(maxn2, diag = T)] = NA
     
     for (projnum in 1:nrow(dfprojns)){
          #projnum = 174
          print(projnum)
          #genepairs <- which(genepairsprojn==projnum, arr.ind = T)
          both = intersect(colnames(treesObj$report)[as.logical(as.numeric(strsplit(dfprojns$projout[projnum],'')[[1]]))],
                           cm)
          nb=length(both)
          #tree1 <- treesObj$trees[[genepairs[1,1]]]
          if(nb>2){
               tree1=unroot(pruneTree(treesObj$master,both))
          }
          #find all the genes that contain all of the species in tree1
          allreport=treesObj$report[,both]
          if(nb>1){
               ss=rowSums(allreport)
          }else{
               ss=allreport
          }
          iiboth=which(ss==length(both))
          ai2=which(maxn2[iiboth, iiboth]==nb, arr.ind = T)
          
          if(length(both)<min.sp){
               #donei[iiboth[ai]]=1
               corout[matrix(iiboth[ai2],nrow=nrow(ai2))]=Inf
               ai2rev <- ai2
               ai2rev[,1] <- ai2[,2]
               ai2rev[,2] <- ai2[,1]
               corout[matrix(iiboth[ai2rev],nrow=nrow(ai2rev))]=1
               pairsdone <- pairsdone+nrow(ai2)
               message(paste0(pairsdone," done out of ",totalpairs))
               next
          }    
          
          if(T){
               
               ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
               
               ii= treesObj$matIndex[ee[, c(2,1)]]
               
               allbranch=treesObj$paths[iiboth,ii]
               if(weighted){
                    allbranchw=weights[iiboth,ii]
               }
               if(scaleForPproj){
                    nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
               }else{
                    nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
               }
               
               iibad=which(allbranch<cutoff)
               #don't scale
               #allbranch=scaleMat(allbranch)
               if(!is.null(transform)){
                    nv=transform(nv)
                    allbranch=transform(allbranch)
               }
               allbranch[iibad]=NA
               
               if(!scale){
                    if(!weighted){
                         proj=residfunc(allbranch, model.matrix(~1+nv))
                         
                    }else{
                         
                         proj=residfunc(allbranch, model.matrix(~1+nv), allbranchw)
                         
                    }
               }else{
                    
                    if(!weighted){
                         proj=residfunc(allbranch, model.matrix(~1+nv))
                    }else{
                         
                         proj=residfunc(allbranch, model.matrix(~1+nv),allbranchw)
                    }
                    
                    proj=scale(proj, center = F)                    
               }
          }
               
          for(m in 1:nrow(ai2)){
               #m = 1
               #print(m)
               k=sort(ai2[m,])[1]
               l=sort(ai2[m,])[2]
               #tmpcor=cor(proj[k,], proj[l,], method = 'p', use = 'pair')
               if(sum(!is.na(proj[k,])&!is.na(proj[l,])) >= min.sp){
                    cres=cor.test(win(proj[k,],3), win(proj[l,],3), method='pearson', exact=F)
                    corout[iiboth[l], iiboth[k]]=cres$est
                    corout[iiboth[k], iiboth[l]]=sum(!is.na(proj[k,])&!is.na(proj[l,]))
               }else{
                    corout[iiboth[l], iiboth[k]]=Inf
                    corout[iiboth[k], iiboth[l]]=sum(!is.na(proj[k,])&!is.na(proj[l,]))
               }
	       #corout[iiboth[k],iiboth[l]] <- 1
               #print(m)
          }
          pairsdone <- pairsdone+nrow(ai2)
          message(paste0(pairsdone," done out of ",totalpairs))
     }
     corout
}
getAllProjections_single <- function(treesrdsfilepath){
     trefile <- gsub('.rds','',tail(strsplit(treesrdsfilepath,'/')[[1]],n =1))
     treeObj <- readRDS(treesrdsfilepath)
     treesreport <- treeObj$report
     rm('treeObj')
     maxn=treesreport%*%t(treesreport)
     ltmaxn <- which(lower.tri(maxn), arr.ind = T)
     rm('maxn')
     print(paste0(nrow(ltmaxn),' gene pairs'))
     print(Sys.time())
     res <- sapply(1:nrow(ltmaxn), function(x){
          if(x %% 1000000 == 0){
               print(x)
               print(Sys.time())
          }
          paste(treesreport[ltmaxn[x,1],]*treesreport[ltmaxn[x,2],], collapse = '')
     })
     print(Sys.time())
     resdf <- data.frame(projout = res, stringsAsFactors = F) 
     rm('res')
     resdf %>% group_by(projout) %>% summarise(nPairs = length(projout)) %>%
          arrange(projout) %>% 
          mutate(nSpecies = sapply(projout,function(y) str_count(y,'1'), USE.NAMES = F))
}
win=function(x,w){
     if(sum(!is.na(x))<4){
          return(x)
     }
     xs=sort(x[!is.na(x)], decreasing = T)
     xmax=xs[w]
     xmin=xs[length(xs)-w+1]
     
     x[x>xmax]=xmax
     x[x<xmin]=xmin
     x
}
correlateERCAllTrees=function(treesObj, cutoff=NULL, transform="sqrt", weighted=T, weights = NULL, useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05){
     
     totalpairs <- choose(treesObj$numTrees,2)
     #      treesObj = dtre
     #      cutoff<-NULL
     #      transform="sqrt"
     #      weighted=T
     #      weights = weights.sqrt
     #      useSpecies=NULL
     #      min.sp=5
     #      scale=T
     #      doOnly=NULL
     #      maxT=NULL
     #      scaleForPproj=F
     #      mean.trim=0.05
     
     
     if(is.null(cutoff)){
          cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
          message(paste("cutoff is set to", cutoff))
     }
     if (weighted){
          if(is.null(weights)){
               message("Weights can not be NULL if weighted = T; computing weights")
               weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
          }
          residfunc=fastLmResidMatWeighted
     }else{
          residfunc=fastLmResidMat
     }
     
     if (is.null(useSpecies)){
          useSpecies=treesObj$masterTree$tip.label
          #mappedEdges=trees$mappedEdges
     }
     if(is.null(maxT)){
          maxT=treesObj$numTrees
     }
     if(transform!="none"){
          transform=match.arg(transform,c("sqrt", "log"))
          transform=get(transform)
     }else{
          transform=NULL
     }
     #cm is the names of species that are included in useSpecies and the master tree
     cm=intersect(treesObj$masterTree$tip.label, useSpecies)
     sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
     if (length(sp.miss) > 0) {
          message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                       collapse = ",")))
          
     }
     
     corout=matrix(NA,treesObj$numTrees,treesObj$numTrees)
     diag(corout) = 1
     
     maxn2=treesObj$report%*%t(treesObj$report)
     maxn2[upper.tri(maxn2, diag = T)] = NA
     
     #      if(is.null(doOnly)){
     #           doOnly=1
     #      }else{
     #           maxT=1
     #      }
     pairsdone <- 0
     print(maxT)
     for (i in 2:maxT){
          for(j in 1:(i-1)){
               #i <- 480
               #j <- 2
               #i <- 160
               #j <- 67
               if(is.na(corout[i,j])){
                    #tree1=treesObj$trees[[i]]
                    #get the common species, prune and unroot
                    both=intersect(intersect(treesObj$trees[[i]]$tip, treesObj$trees[[j]]$tip), cm)
                    nb=length(both)
                    if(nb>2){
                         tree1=unroot(pruneTree(treesObj$master,both))
                    }
                    #find all the genes that contain all of the species in tree1
                    allreport=treesObj$report[,both]
                    if(nb>1){
                         ss=rowSums(allreport)
                    }else{
                         ss=allreport
                    }
                    iiboth=which(ss==length(both))
                    ai2=which(maxn2[iiboth, iiboth]==nb, arr.ind = T)
                    
                    if(length(both)<min.sp){
                         corout[matrix(iiboth[ai2],nrow=nrow(ai2))]=Inf
                         ai2rev <- ai2
                         ai2rev[,1] <- ai2[,2]
                         ai2rev[,2] <- ai2[,1]
                         corout[matrix(iiboth[ai2rev],nrow=nrow(ai2rev))]=1
                         pairsdone <- pairsdone+nrow(ai2)
                         message(paste0(pairsdone," done out of ",totalpairs))
                         next
                    }    
                    
                    
                    if(T){
                         
                         
                         ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
                         
                         ii= treesObj$matIndex[ee[, c(2,1)]]
                         
                         allbranch=treesObj$paths[iiboth,ii]
                         if(weighted){
                              allbranchw=weights[iiboth,ii]
                         }
                         if(scaleForPproj){
                              nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
                         }else{
                              nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
                         }
                         
                         iibad=which(allbranch<cutoff)
                         #don't scale
                         #allbranch=scaleMat(allbranch)
                         if(!is.null(transform)){
                              nv=transform(nv)
                              allbranch=transform(allbranch)
                         }
                         allbranch[iibad]=NA
                         
                         if(!scale){
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                                   
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv), allbranchw)
                                   
                              }
                         }else{
                              
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv),allbranchw)
                              }
                              
                              proj=scale(proj, center = F)
                              
                         }
                         
                         
                         #we have the projection
                         
                         for(m in 1:nrow(ai2)){
                              #m = 1
                              #print(m)
                              k=sort(ai2[m,])[1]
                              l=sort(ai2[m,])[2]
                              #tmpcor=cor(proj[k,], proj[l,], method = 'p', use = 'pair')
                              if(sum(!is.na(proj[k,])&!is.na(proj[l,])) >= min.sp){
                                   cres=cor.test(win(proj[k,],3), win(proj[l,],3), method='pearson', exact=F)
                                   if(!is.na(cres$est)){
                                        corout[iiboth[l], iiboth[k]]=cres$est
                                   }else{
                                        corout[iiboth[l], iiboth[k]]=Inf
                                   }
                              }else{
                                   corout[iiboth[l], iiboth[k]]=Inf
                              }
                              corout[iiboth[k], iiboth[l]]=sum(!is.na(proj[k,])&!is.na(proj[l,]))
                         }
                         pairsdone <- pairsdone+nrow(ai2)
                         message(paste0(pairsdone," done out of ",totalpairs))
                    }
                    
               }
               
          }
     }
     rownames(corout) <- names(treesObj$trees)
     colnames(corout) <- names(treesObj$trees)
     corout[is.infinite(corout)] <- NA
     corout
}

