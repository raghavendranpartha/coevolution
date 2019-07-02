library(RERconverge)
require(dplyr)
require(data.table)
hsap.genes <- readLines('../data/hsap.genes')
hsap.vert.genes <- readLines('../data/hsap.vert.genes')
hsap.dmel.genes <- readLines('../data/hsap.dmel.genes')
hsap.cele.genes <- readLines('../data/hsap.cele.genes')
hsap.scer.genes <- readLines('../data/hsap.scer.genes')

genepresence <- (hsap.genes %in% hsap.genes)+(hsap.genes %in% hsap.vert.genes)+(hsap.genes %in% hsap.dmel.genes)+(hsap.genes %in% hsap.scer.genes)+(hsap.genes %in% hsap.cele.genes)
names(genepresence) <- hsap.genes
genepresenceeach <- rbind((hsap.genes %in% hsap.genes),(hsap.genes %in% hsap.vert.genes),(hsap.genes %in% hsap.dmel.genes),(hsap.genes %in% hsap.cele.genes),(hsap.genes %in% hsap.scer.genes))
rownames(genepresenceeach) = c('mamm','vert','dmel','cele','scer')
colnames(genepresenceeach) = hsap.genes
genepresenceeach[,1:5]

genesinatleast2 <- hsap.genes[which(genepresence>=2)]
nonsinglegenepairs <- as.data.frame(t(combn(genesinatleast2,2)), stringsAsFactors = F)
nonsinglegenepairs$gene1c = genepresence[nonsinglegenepairs$V1]
nonsinglegenepairs$gene2c = genepresence[nonsinglegenepairs$V2]
nonsinglegenepairs$pairc <- with(nonsinglegenepairs, gene1c+gene2c)
nonsinglegenepairs <- arrange(nonsinglegenepairs, desc(pairc))

saveRDS(list(genepresenceeach = genepresenceeach,
             nonsinglegenepairs = nonsinglegenepairs), file = '../data/fulldatasetintegration.rds')
#head(nonsinglegenepairs)

#calc common pairs - h2p

# fulldata <- readRDS('../data/fulldatasetintegration.rds')
# nrs <- nrow(fulldata$nonsinglegenepairs)
# id <- 905
# # rowstodost <- 100000*(id-1)+1
# rowstodoen <- 100000*id
# if(rowstodoen > nrs){
#      rowstodoen = nrs
# }
# #rowstodost <- 1
# #rowstodoen <- 10000
# pairc = sapply(rowstodost:rowstodoen, function(x){
#       if(x%%1000 == 0){print(x)}
#       with(fulldata,sum(genepresenceeach[,nonsinglegenepairs$V1[x]]*genepresenceeach[,nonsinglegenepairs$V2[x]]))
# })
# saveRDS(pairc, file =paste0('../data/paircs/pairc.',id,'.rds'))

fulldata <- readRDS('../data/fulldatasetintegration.rds')

paircs <- c()
for(ii in c(1:906)){
     paircs <- c(paircs, readRDS(paste0('../data/paircs/pairc.',ii,'.rds')))
}

fulldata$nonsinglegenepairs$pairc <- paircs
rm('paircs')
table(fulldata$nonsinglegenepairs$pairc)

saveRDS(fulldata, file = '../data/fulldatasetintegration_correct.rds')

fulldata <- readRDS('../data/fulldatasetintegration_correct.rds')
# fulldata$genepresenceeach[,'KCNRG']
# fulldata$genepresenceeach[4,'KCNRG'] <- FALSE
# head(filter(fulldata$nonsinglegenepairs,V1 == 'KCNRG'))
# head(filter(fulldata$nonsinglegenepairs,V2 == 'KCNRG'))



#testp <- sample_n(filter(fulldata$nonsinglegenepairs, pairc == 5), 100)
testp <-with(fulldata,rbind(sample_n(filter(nonsinglegenepairs, pairc == 5), 100),
      sample_n(filter(nonsinglegenepairs, pairc == 4), 100),
      sample_n(filter(nonsinglegenepairs, pairc == 3), 100),
      sample_n(filter(nonsinglegenepairs, pairc == 2), 100),
      sample_n(filter(nonsinglegenepairs, pairc == 1), 100)))

sum(nonsinglegenepairs$pairc >= 3)
sum(nonsinglegenepairs$pairc == 2)
sum(nonsinglegenepairs$pairc == 1)
sum(nonsinglegenepairs$pairc == 0)


getercasvec <- function(mat, gen){
     #mat <- merc
     genv <- match(gen,rownames(mat))
     genord <- gen[order(genv)]
     out <- mat[genord,genord][lower.tri(mat[genord,genord], diag = F)]
     #out[is.infinite(out)] <- NA
     out
}

#integration

trs <- c('~/Documents/erc/data/saves/mamm63nt.trees.rds', '~/Documents/erc/data/saves/vert39.trees.rds', '~/Documents/erc/data/saves/alldroso22.tre.rds' ,'~/Documents/erc/data/saves/worms17.trees.rds', '~/Documents/erc/data/saves/all_Scer.tre.rds')
ercs <- c('~/Documents/erc/data/saves/mamm63.erc.rds', '~/Documents/erc/data/saves/vert39.erc.rds', '~/Documents/erc/data/saves/alldroso22.tre.erc.rds', '~/Documents/erc/data/saves/worms17.erc.rds','~/Documents/erc/data/saves/all_Scer.erc.rds')
betas <- c('../data/ercbetafit/mamm63.erc.betafits.rds', '../data/ercbetafit/vert39.erc.betafits.rds', '../data/ercbetafit/alldroso22.tre.erc.betafits.rds', '../data/ercbetafit/worms17.erc.betafits.rds','../data/ercbetafit/all_Scer.erc.betafits.rds')
mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')


#for(ii in c(1:5)){
ii <- 3 #mamm,vert,dmel,cele,scer
tre <- readRDS(trs[ii])
erc <- readRDS(ercs[ii])
betafit <- readRDS(betas[ii])
betaspranges <- cbind(as.numeric(sapply(strsplit(names(betafit),'_'),'[[',1)),
                      as.numeric(sapply(strsplit(names(betafit),'_'),'[[',2)))
erc[is.infinite(erc)] <- NA
maxn <- tre$report%*%t(tre$report)
maxn[upper.tri(maxn, diag = T)] = NA
vals <- t(sapply(1:nrow(testp), function(x){
     #x <- 1
     g1 <- testp$V1[x]
     g2 <- testp$V2[x]
     g1 <- "CISD3"
     g2 <- "KCNRG"
     if(ii < 3){
          ercv <- ifelse(fulldata$genepresenceeach[ii,g1]&fulldata$genepresenceeach[ii,g2],
                        getercasvec(erc,c(g1,g2)),NA)
          nspv <- ifelse(fulldata$genepresenceeach[ii,g1]&fulldata$genepresenceeach[ii,g2],
                         getercasvec(maxn,c(g1,g2)),NA)
     }else{
          ercv <- ifelse(fulldata$genepresenceeach[ii,g1]&fulldata$genepresenceeach[ii,g2],
                        getercasvec(erc,c(mapvecs[[ii-2]][g1],mapvecs[[ii-2]][g2])),NA)
          nspv <- ifelse(fulldata$genepresenceeach[ii,g1]&fulldata$genepresenceeach[ii,g2],
                         getercasvec(maxn,c(mapvecs[[ii-2]][g1],mapvecs[[ii-2]][g2])),NA)
     }
     if(!is.na(ercv)){
          if(mapvecs[[ii-2]][g1]==mapvecs[[ii-2]][g2]){
               ercv <- NA
               nspv <- NA
          }else{
               whichbetafit <- which((betaspranges <= nspv)[,1] & (betaspranges > nspv)[,2])
          }
     }
     #whichbetafit <- which((betaspranges <= nspv)[,1] & (betaspranges > nspv)[,2])
     if(!is.na(ercv)){
          pv1 <- pbeta((ercv+1)/2,betafit[[whichbetafit]][1],betafit[[whichbetafit]][2],lower.tail = F)
          pv2 <- pbeta((ercv+1)/2,betafit[[whichbetafit]][1],betafit[[whichbetafit]][2],lower.tail = T)
     }else{
          pv1 <- NA
          pv2 <- NA
     }
     c(ercv,nspv,-log(pv1,10),-log(pv2,10),max(c(-log(pv1,10),-log(pv2,10))))
}))
colnames(vals) <- paste0(c('erc','nsp','nlogpv1','nlogpv2','nlogpvbest'),ii)
#}
testps <- cbind(testp,vals)

fd <- readRDS('../data/integration/int1_5.rds')
for(ii in c(2:5)){
     fd <- rbind(fd, readRDS(paste0('../data/integration/int',ii,'_5.rds')))
}

mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')

