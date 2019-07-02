require(dplyr)
require(data.table)

library(topGO)
library(GOSemSim)
library(org.Hs.eg.db)

xx <- as.list(org.Hs.egALIAS2EG)
toppairs <- fread('../data/toperc/sumlogpftercgt15.bestbitscorepclt10.tsv', header = T)
toppairs$genesimbp <- NA
nrows <- 100
toppairs$genesimbp[1:nrows] <- getgenesim.bp(toppairs, nrows)
toppairs$genesimbp[is.infinite(toppairs$genesimbp)] <- NA
quantile(toppairs$genesimbp, na.rm = T)
getgenesim.bp <- function(fd,nrows){
     genesim.bp <- sapply(1:nrows, function(ii){
          print(ii)
          g1 <- fd$g1[ii]
          g2 <- fd$g2[ii]
          eg1 <- ifelse(g1 %in% names(xx), xx[[g1]], NA) 
          eg2 <- ifelse(g2 %in% names(xx), xx[[g2]], NA)
          if(!is.na(eg1) & !is.na(eg2)){
               combos <- expand.grid(eg1,eg2)
               gs <-max(sapply(1:nrow(combos), function(x){
                    x <- 1
                    gslist <- geneSim(combos$Var1[x],combos$Var2[x], ont = 'BP')
                    if(!is.na(gslist)){
                         gslist$geneSim
                    }else{
                         NA
                    }
               }), na.rm = T)
          }
          gs
     })
}