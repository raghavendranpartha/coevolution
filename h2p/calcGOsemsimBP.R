args <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(data.table)
library(GOSemSim)
library(org.Hs.eg.db)

inf <- args[1]
njobs <- as.numeric(args[2])
split <- as.numeric(args[3])
#inf <- '../data/toperc/int.sumnlogpv1gteq2paircgteq3bitscorepclt10.tsv'
#njobs <- 1000
#split <- 2

fn <- gsub('.tsv','',tail(strsplit(inf,'/')[[1]],1))
int.min2.pairc3.dup <- fread(inf, header = T)
xx <- as.list(org.Hs.egALIAS2EG)
d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)

vsplit <- function(v, n) {
     l = length(v)
     r = l/n
     return(lapply(1:n, function(i) {
          s = max(1, round(r*(i-1))+1)
          e = min(l, round(r*i))
          return(v[s:e])
     }))
}

out <- vsplit(1:nrow(int.min2.pairc3.dup), njobs)
start = out[[split]][1]
end = tail(out[[split]],1)

getgenesim.bp <- function(fd,rowstart, rowend){
     genesim.bp <- sapply(rowstart:rowend, function(ii){
          print(ii)
          g1 <- fd$V1[ii]
          g2 <- fd$V2[ii]
          eg1 <- ifelse(g1 %in% names(xx), xx[[g1]], NA) 
          eg2 <- ifelse(g2 %in% names(xx), xx[[g2]], NA)
          if(!is.na(eg1) & !is.na(eg2)){
               combos <- expand.grid(eg1,eg2)
               gs <-max(sapply(1:nrow(combos), function(x){
                    x <- 1
                    gslist <- geneSim(combos$Var1[x],combos$Var2[x], semData = d)
                    if(!is.na(gslist)){
                         gslist$geneSim
                    }else{
                         NA
                    }
               }), na.rm = T)
          }
          gs
     })
     genesim.bp[is.infinite(genesim.bp)] <- NA
     genesim.bp
}

genesimbp <- getgenesim.bp(int.min2.pairc3.dup, start, end)
outdir <- '../data/geneSimBP/'
dir.create(outdir, showWarnings = F)
saveRDS(genesimbp, paste0(outdir,fn,'.',split,'.rds'))
