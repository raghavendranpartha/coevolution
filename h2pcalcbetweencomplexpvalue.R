args <- commandArgs(trailingOnly = T)

require(dplyr)
require(data.table)
require(tidyr)
source('../code/funcsERCfromProjections.R')

cplxf <- args[1]
cplxpairf <- args[2]
ii <- args[3]

outf <- gsub('.rds','',tail(strsplit(cplxf,'/',fixed = T)[[1]],1))
resdir <- paste0('../data/randoms/betweencomplex/',outf)

corum.cplxsq <- fread(cplxpairf)
corum.cplxsq$pair <- with(corum.cplxsq, paste0(complexind1,'_',complexind2))
sumpv <- npv <- rep(0,length(corum.cplxsq$pair))
names(sumpv) <- corum.cplxsq$pair
names(npv) <- corum.cplxsq$pair

dff <- fread(paste0(resdir,'/runs',ii,'.tsv'), header = T) %>% mutate(pair = paste0(complexind1,'_',complexind2))
dffg <- spread(select(dff,c('pair','xv','dataset')), key = dataset, value=xv)
dffgmat <- as.matrix(dffg[,c(2:11)])
rownames(dffgmat) <- dffg$pair
dffgmatu <- dffgmat[corum.cplxsq$pair,]
rm('dff');rm('dffg');rm('dffgmat')
sumpv <- sumpv+rowSums(dffgmatu>=corum.cplxsq$meaninterc,na.rm = T)
npv <- npv+rowSums(!is.na(dffgmatu))

saveRDS(list(sumpv=sumpv,npv=npv),file=paste0(resdir,'/runs',ii,'.rds'))

sumpv <- readRDS('../data/randoms/betweencomplex/corumcplxs.cleaned/runs1.rds')$sumpv
npv <- readRDS('../data/randoms/betweencomplex/corumcplxs.cleaned/runs1.rds')$npv
for(ii in 2:1000){
     df <- readRDS(paste0('../data/randoms/betweencomplex/corumcplxs.cleaned/runs',ii,'.rds'))
     if(!identical(names(df$sumpv),names(sumpv))){
          print(ii)
     }else{
          sumpv <- sumpv+df$sumpv
          npv <- npv+df$npv
     }
}
corum.cplxsq$pv <- sumpv/npv
corum.cplxsq$qv <- p.adjust(corum.cplxsq$pv,method = 'BH')
corum.cplxsq <- arrange(corum.cplxsq,qv)
sum(corum.cplxsq$qv < 0.05)
View(filter(corum.cplxsq,qv<0.05))
length(unique(c(filter(corum.cplxsq,qv<0.05)$complexind1,
                filter(corum.cplxsq,qv<0.05)$complexind2)))
overrepcomplexes <- c(filter(corum.cplxsq,qv<0.05)$complex1,
                  filter(corum.cplxsq,qv<0.05)$complex2)
head(sort(table(overrepcomplexes), decreasing = T),5)
tail(sort(table(overrepcomplexes), decreasing = T))
write.table(corum.cplxsq,file = '../data/public/betweencomplex/corum.betweencomplex.withpqvalues.tsv',
            quote = F, sep = '\t', row.names = F)