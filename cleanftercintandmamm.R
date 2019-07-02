badblastpairs <- fread('../data/allgeneshg19.pairbitscore.gteq0p1.tsv')

sumfterc.raw <- readRDS('../data/RAWallsum.fterc.rds')#from the server output

sum(sumfterc.raw == 0)
sumfterc.raw[sumfterc.raw == 0] <- NA
diag(sumfterc.raw) <- NA
sum(is.na(sumfterc.raw))# = 6428897, ~6M

for(ii in 1:nrow(badblastpairs)){
     if(ii %% 100 == 0){
          print(ii)
     }
     sumfterc.raw[badblastpairs$g1[ii],badblastpairs$g2[ii]] <- NA
     sumfterc.raw[badblastpairs$g2[ii],badblastpairs$g1[ii]] <- NA
}

saveRDS(sumfterc.raw, '../data/allsum.fterc.rds')

merc.raw <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds')#from the server output
merc.raw[upper.tri(merc.raw)] <- 0
merc.raw <- merc.raw+t(merc.raw)
diag(merc.raw) <- NA
sum(is.na(merc.raw))# ~ 10M

for(ii in 1:nrow(badblastpairs)){
     if(ii %% 100 == 0){
          print(ii)
     }
     merc.raw[badblastpairs$g1[ii],badblastpairs$g2[ii]] <- NA
     merc.raw[badblastpairs$g2[ii],badblastpairs$g1[ii]] <- NA
}

saveRDS(merc.raw, '~/Documents/erc/data/saves/mamm63nt.trees.fterc.symm.rds')

