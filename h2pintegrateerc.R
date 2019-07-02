args <- commandArgs(trailingOnly=TRUE)

trs <- c('/zfs1/nclark/rap119/erc/data/trees/mamm63nt.trees.rds', 
         '/zfs1/nclark/rap119/erc/data/trees/vert39.trees.rds', 
         '/zfs1/nclark/rap119/erc/data/trees/alldroso22.tre.rds',
         '/zfs1/nclark/rap119/erc/data/trees/worms17.trees.rds', 
         '/zfs1/nclark/rap119/erc/data/trees/all_Scer.tre.rds')
ercs <- c('/zfs1/nclark/rap119/erc/data/saves/ercmat/mamm63nt.trees.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/ercmat/vert39.trees.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/ercmat/alldroso22.tre.erc.rds', 
          '/zfs1/nclark/rap119/erc/data/saves/ercmat/worms17.trees.erc.rds',
          '/zfs1/nclark/rap119/erc/data/saves/ercmat/all_Scer.tre.erc.rds')
betas <- c('/zfs1/nclark/rap119/erc/data/saves/ercbetafit/mamm63nt.trees.erc.betafits.rds', 
           '/zfs1/nclark/rap119/erc/data/saves/ercbetafit/vert39.trees.erc.betafits.rds', 
           '/zfs1/nclark/rap119/erc/data/saves/ercbetafit/alldroso22.tre.erc.betafits.rds', 
           '/zfs1/nclark/rap119/erc/data/saves/ercbetafit/worms17.trees.erc.betafits.rds',
           '/zfs1/nclark/rap119/erc/data/saves/ercbetafit/all_Scer.tre.erc.betafits.rds')
mapvecs <- readRDS('/zfs1/nclark/rap119/coevolution/data/flywormyeastmap.vec.rds')
minspv <- c(15,10,10,10,10)

fulldata <- readRDS(file = '/zfs1/nclark/rap119/coevolution/data/fulldatasetintegration_correct.rds')
id <- as.numeric(args[1])
rowstodost <- 100000*(id-1)+1
rowstodoen <- 100000*id
nrs <- nrow(fulldata$nonsinglegenepairs)
if(rowstodoen > nrs){
      rowstodoen <- nrs
}

genepresenceeach = fulldata$genepresenceeach
testp <- fulldata$nonsinglegenepairs[rowstodost:rowstodoen,]

rm('fulldata')

ii <- as.numeric(args[2])
#ii <- 2

getercasvec <- function(mat, gen){
     #mat <- merc
     genv <- match(gen,rownames(mat))
     genord <- gen[order(genv)]
     out <- mat[genord,genord][lower.tri(mat[genord,genord], diag = F)]
     #out[is.infinite(out)] <- NA
     out
}
getnspasvec <- function(mat, gen){
     #mat <- merc
     genv <- match(gen,rownames(mat))
     genord <- gen[order(genv, decreasing = T)]
     out <- mat[genord,genord][lower.tri(mat[genord,genord], diag = F)]
     #out[is.infinite(out)] <- NA
     out
}


tre <- readRDS(trs[ii])
erc <- readRDS(ercs[ii])
betafit <- readRDS(betas[ii])
betaspranges <- cbind(as.numeric(sapply(strsplit(names(betafit),'_'),'[[',1)),
                      as.numeric(sapply(strsplit(names(betafit),'_'),'[[',2)))
vals <- t(sapply(1:nrow(testp), function(x){
     #x <- 1
     g1 <- testp$V1[x]
     g2 <- testp$V2[x]
     print(g1)
     print(g2)
     if(ii < 3){#this is for mammalian or vertebrate beta pvalue calculations
          ercv <- ifelse(genepresenceeach[ii,g1]&genepresenceeach[ii,g2],
                        getercasvec(erc,c(g1,g2)),NA)
          nspv <- ifelse(genepresenceeach[ii,g1]&genepresenceeach[ii,g2],
                         getnspasvec(erc,c(g1,g2)),NA)
     }else{
          ercv <- ifelse(genepresenceeach[ii,g1]&genepresenceeach[ii,g2],
                        getercasvec(erc,c(mapvecs[[ii-2]][g1],mapvecs[[ii-2]][g2])),NA)
          nspv <- ifelse(genepresenceeach[ii,g1]&genepresenceeach[ii,g2],
                         getnspasvec(erc,c(mapvecs[[ii-2]][g1],mapvecs[[ii-2]][g2])),NA)
     }
     #nspv <- 22
     if(!is.na(ercv)){
          if(ii>2){
               if(mapvecs[[ii-2]][g1]==mapvecs[[ii-2]][g2]){
                    ercv <- NA
                    nspv <- NA
               }else{
                    whichbetafit <- which((betaspranges <= nspv)[,1] & (betaspranges > nspv)[,2])
               }
          }else{
               whichbetafit <- which((betaspranges <= nspv)[,1] & (betaspranges > nspv)[,2])
          }
     }
     if(!is.na(nspv)){
          if(nspv < minspv[ii]){
               ercv <- NA
               nspv <- NA
          }
     }
     #whichbetafit <- which((betaspranges <= nspv)[,1] & (betaspranges > nspv)[,2])
     if(!is.na(ercv)){
          pv1 <- pbeta((ercv+1.01)/2.02,betafit[[whichbetafit]][1],betafit[[whichbetafit]][2],lower.tail = F)
          pv2 <- pbeta((ercv+1.01)/2.02,betafit[[whichbetafit]][1],betafit[[whichbetafit]][2],lower.tail = T)
     }else{
          pv1 <- NA
          pv2 <- NA
     }
     c(ercv,nspv,-log(pv1,10),-log(pv2,10),max(c(-log(pv1,10),-log(pv2,10))))
}))
colnames(vals) <- paste0(c('erc','nsp','nlogpv1','nlogpv2','nlogpvbest'),ii)

testps <- cbind(testp,vals)
saveRDS(testps, file = paste0('../data/integration/int',ii,'/int',id,'_',ii,'.rds'))
