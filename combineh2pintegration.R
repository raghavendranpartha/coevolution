require(dplyr)
#require(data.table)

#fd <- readRDS('../data/fulldatasetintegration_correct.rds')$nonsinglegenepairs
#fd[,6:30] <- matrix(NA, nrow = nrow(fd),ncol = 25)

#jj <- 1
#fd <- readRDS(paste0('../data/integration/int',jj,'/int',jj,'_1.rds'))
# thr <- 3
# fd <- readRDS(file = paste0('../data/integration/intall/int',1,'.rds')) %>%
#      #filter(sumnlogpvbest >= 3)
#      filter(sumnlogpv1 >= thr)
# for(ii in c(2:906)){
#      print(ii)
#      fd <- rbind(fd,readRDS(file = paste0('../data/integration/intall/int',ii,'.rds')) %>%
#           filter(sumnlogpvbest >= thr))
# }
# write.table(fd, file = '../data/toperc/int.sumnlogpvbest.gteq.10.tsv',
#             quote = F, sep = '\t',row.names = F)

x<-1
fd <- readRDS(paste0('../data/integration/intall/int',x,'.rds'))
thr <- 2
alldf <- lapply(1:906, function(x){
     print(x)
     readRDS(paste0('../data/integration/intall/int',x,'.rds')) %>% 
          filter(sumnlogpv1 >= thr)
})
fd <- do.call(rbind, alldf)

write.table(fd, file = paste0('../data/toperc/int.sumnlogpv1.gteq.',thr,'.tsv'),
            quote = F, sep = '\t',row.names = F)

thr <- 0.2
alldf <- lapply(1:906, function(x){
     print(x)
     readRDS(paste0('../data/integration/intall/int',x,'.rds')) %>% 
          filter(sumnlogpv1 < thr)
})
fd <- do.call(rbind, alldf)

write.table(fd, file = paste0('../data/toperc/int.sumnlogpv1.lt.',thr,'.tsv'),
            quote = F, sep = '\t',row.names = F)



     
# fd <- readRDS('../data/int99_5.rds') %>%
#      na.omit()
# View(fd[1:1000,])

#vis egs
# library(RERconverge)
# source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
# source('~/Documents/RERc/RERconverge/R/RcppExports.R')
# source('~/Documents/RERc/RERconverge/R/projection_coevo.R')
# source('~/Documents/erc/code/funcsERCfromProjections.R')
# require(gridExtra)
# gns <- c('TSN','MYL12A')
# mtre <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
# mwts <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.weights.rds')
# mrers <- readRDS('~/Documents/rermethods/data/mamm63nt.trees.scaledrers.sqrt.wt.rds')
# correlateERCGeneListTrees(mtre, gns, weighted = T, weights = mwts)
# grid.arrange(plotRers(mrers,'TSN'),plotRers(mrers,'MYL12A'))



for(ii in c(1:906)){
     print(ii)
     #ii <- 1
     #jj <- 1
     for(jj in c(1:5)){
          assign(paste0('fd',jj), readRDS(paste0('../data/integration/int',jj,'/int',ii,'_',jj,'.rds')))
     }
     fd <- cbind(fd1,fd2[,6:10],fd3[,6:10],fd4[,6:10],fd5[,6:10]) 
     fd$sumnlogpv1 <- rowSums(fd[,c('nlogpv11','nlogpv12','nlogpv13','nlogpv14','nlogpv15')], na.rm = T)
     fd$sumnlogpvbest <- rowSums(fd[,c('nlogpvbest1','nlogpvbest2','nlogpvbest3','nlogpvbest4','nlogpvbest5')], na.rm = T)
     fd <- select(fd,c(1:5,31,32,6:30)) %>%
          arrange(desc(sumnlogpvbest))
     saveRDS(fd, file = paste0('../data/integration/intall/int',ii,'.rds'))
}

#View(fd[1:5,])
