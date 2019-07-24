require(dplyr)
require(data.table)
require(ggplot2)

for(onttype in c('CC','MF','BP')){
#for(onttype in c('CC','MF')){
     #onttype <- 'BP'
     assign(paste0('top100k',onttype),
            readRDS(paste0('../data/geneSim/',onttype,'/top.sumfterc5.bestbitscorepclt10.tsv.rds')))
     assign(paste0('bgd500k',onttype),
            readRDS(paste0('../data/geneSim/',onttype,'/controlndatasetsmatched.top.sumfterc5.bestbitscorepclt10.tsv.rds')))
}


toppairs <- fread('../data/toperc/top.sumfterc5.bestbitscorepclt10.tsv', header = T) %>%
     mutate(pair = as.character(pair)) %>% filter(pair %in% names(top100kCC))
controltoppairs <- fread('../data/toperc/controlndatasetsmatched.top.sumfterc5.bestbitscorepclt10.tsv') %>%
     mutate(pair = as.character(pair)) %>%
     filter(pair %in% names(bgd500kCC))

toppairs$mfsemsim <- top100kMF[toppairs$pair]
controltoppairs$mfsemsim <- bgd500kMF[controltoppairs$pair]
toppairs$bpsemsim <- top100kBP[toppairs$pair]
controltoppairs$bpsemsim <- bgd500kBP[controltoppairs$pair]
toppairs$ccsemsim <- top100kCC[toppairs$pair]
controltoppairs$ccsemsim <- bgd500kCC[controltoppairs$pair]

controltoppairs.control <- filter(controltoppairs, interc > -2, interc < 2)
controltoppairs.control <- controltoppairs

ns <- c(1000,10000,100000)
for(nn in ns){
     mfp <- wilcox.test(toppairs[1:nn,]$mfsemsim, controltoppairs.control$mfsemsim, alternative = 'g')$p.value
     ccp <- wilcox.test(toppairs[1:nn,]$ccsemsim, controltoppairs.control$ccsemsim, alternative = 'g')$p.value
     bpp <- wilcox.test(toppairs[1:nn,]$bpsemsim, controltoppairs.control$bpsemsim, alternative = 'g')$p.value
     print(paste0(nn,':',formatC(bpp, 3),',',formatC(ccp, 3),',',format(mfp, 3)))
}

bpx <- qqplot(controltoppairs.control$ccsemsim, toppairs[1:10000,]$bpsemsim, col = 'red',plot = F)
bpx2 <- qqplot(controltoppairs.control$ccsemsim, toppairs[1:1000,]$bpsemsim, col = 'red',plot = F)
mfx <- qqplot(controltoppairs.control$mfsemsim, toppairs[1:10000,]$mfsemsim, col = 'red',plot = F)
mfx2 <- qqplot(controltoppairs.control$mfsemsim, toppairs[1:1000,]$mfsemsim, col = 'red',plot = F)
ccx <- qqplot(controltoppairs.control$ccsemsim, toppairs[1:10000,]$ccsemsim, col = 'red',plot = F)
ccx2 <- qqplot(controltoppairs.control$ccsemsim, toppairs[1:1000,]$ccsemsim, col = 'red',plot = F)

quantile(toppairs$ccsemsim, probs = seq(0,1,0.1), na.rm = T)
quantile(controltoppairs.control$ccsemsim, probs = seq(0,1,0.1), na.rm = T)

hist(controltoppairs.control$bpsemsim, breaks = 25, probability = T, col = rgb(0,0,0,0.25))
hist(toppairs$bpsemsim, breaks = 25, probability = T, add = T, col = rgb(1,0,0,0.25))

pdf('../figures/gosemsim.pdf',height = 4, width = 8,family='FreeSans')
par(mfrow=c(1,3))
qqplot(controltoppairs.control$bpsemsim, toppairs$bpsemsim, xlab = 'Control pairs', ylab = 'Top pairs', main = '', type = 'l',
       lwd = 2,cex.lab=1.25)
title('A. Biological Pathway', adj = 0, font.main = 1, cex.main = 1.5)
lines(bpx[['x']],bpx[['y']], col = 'red',lwd = 2)
lines(bpx2[['x']],bpx2[['y']], col = 'blue',lwd = 2)
legend('bottomright',c('Top 1000','Top 10,000','Top 100,000'), col = c('blue','red','black'),bty = 'n', pch = 19, cex = 1.25,
       title = 'Dataset')
abline(a=0,b=1,lty=2,lwd=2)
qqplot(controltoppairs.control$ccsemsim, toppairs[1:100000,]$ccsemsim, xlab = 'Control pairs', ylab = 'Top pairs', main = '',
       type = 'l', lwd = 2, cex.lab = 1.25)
title('B. Cellular Component', adj = 0, font.main = 1, cex.main = 1.5)
lines(ccx[['x']],ccx[['y']], col = 'red', lwd = 2)
lines(ccx2[['x']],ccx2[['y']], col = 'blue', lwd = 2)
legend('bottomright',c('Top 1000','Top 10,000','Top 100,000'), col = c('blue','red','black'),bty = 'n', pch = 19, cex = 1.25,
       title = 'Dataset')
abline(a=0,b=1, lty = 2, lwd = 2)
qqplot(controltoppairs.control$mfsemsim, toppairs$mfsemsim, xlab = 'Control pairs', ylab = 'Top pairs', main = '', type = 'l',
       lwd = 2, cex.lab = 1.25)
title('C. Molecular Function', adj = 0, font.main = 1, cex.main = 1.5)
lines(mfx[['x']],mfx[['y']], col = 'red', lwd = 2)
lines(mfx2[['x']],mfx2[['y']], col = 'blue', lwd = 2)
legend('bottomright',c('Top 1000','Top 10,000','Top 100,000'), col = c('blue','red','black'),bty = 'n', pch = 19, cex = 1.25,
       title = 'Dataset')
abline(a=0,b=1, lty = 2, lwd = 2)
dev.off()




plotbox2(controls = onembgd.control$mfsemsim,
         toppairs$interc, toppairs$mfsemsim, ylimi = c(0,1.1))
plotbox2(controls = onembgd.control$ccsemsim,
         toppairs$interc, toppairs$ccsemsim, ylimi = c(0,0.4))
plotbox2(controls = onembgd.control$bpsemsim,
         toppairs$interc, toppairs$bpsemsim, ylimi = c(0,1.1))


onembgd$interclv <- cut(onembgd$interc, breaks = c(-6,-2,2,24))
g <- ggplot(filter(onembgd, pair!=0), aes(x=interclv,y=mfsemsim,fill=interclv))+
     geom_boxplot()
print(g)
quantile(bgd, probs = seq(0,1,0.1), na.rm = T)
quantile(top1k, probs = seq(0,1,0.1), na.rm = T)
top1k <- readRDS('../data/geneSimBP/top1k.genesim.bp.rds')

hist(top1k, breaks = 20, col = rgb(1,0,0,alpha = 0.25), probability = T)
hist(bgd, breaks = 20, col = rgb(0,0,0,alpha = 0.25), add = T, probability = T)

plotbox2(bgd, xx = toppairs[1:10000,]$interc, yy = toppairs[1:10000,]$mfsemsim, ylimi = c(0,1))


sum(is.na(bgd))
sum(is.na(top1k))

badcodes <- c('IBA','IBD','IKR','IRD','ISS','ISO','ISA','ISM','IGC','RCA')
table(gogaf$V7)[badcodes]

gogaf2 <- fread('../data/public/GOgaf/RAWgoa_human.gaf', skip = 31, header = F) %>% dplyr::select(c(2,3,5,7,10)) %>%
     filter(!V7 %in% badcodes)
gogaf2.bygene <- group_by(gogaf2, V3) %>% summarise(terms = paste(unique(V5), collapse = ","), nterms = length(unique(V5)))
gogaf2.byterm <- group_by(gogaf2, V5) %>% summarise(genes = paste(unique(V3), collapse = ","), ngenes = length(unique(V3)))
gogaf2.filt <- filter(gogaf2, !V5 %in% filter(gogaf2.byterm,ngenes>500)$V5)
write.table(gogaf2.filt, file = '../data/public/GOgaf/goa_human.onlyexperimental.gaf',
            quote = F, sep = '\t', row.names = F)


gogaf <- fread('../data/public/GOgaf/RAWgoa_human.gaf', skip = 31, header = F) %>% dplyr::select(c(2,3,5,10))
gogaf.bygene <- group_by(gogaf, V3) %>% summarise(terms = paste(unique(V5), collapse = ","), nterms = length(unique(V5)))
gogaf.byterm <- group_by(gogaf, V5) %>% summarise(genes = paste(unique(V3), collapse = ","), ngenes = length(unique(V3)))
gogaf.filt <- filter(gogaf, !V5 %in% filter(gogaf.byterm,ngenes>500)$V5)
write.table(gogaf.filt, file = '../data/public/GOgaf/goa_human.gaf',
            quote = F, sep = '\t', row.names = F)
#gogaf.filt.byterm <- group_by(gogaf.filt, V5) %>% summarise(genes = paste(unique(V3), collapse = ","), ngenes = length(unique(V3)))
gogaf.unip.sym <- dplyr::select(gogaf, c(1,2)) %>% unique()
sum(gogaf.unip.sym$V2 == gogaf.unip.sym$V3)
gogaf.unip.sym.u <- gogaf.unip.sym[!(gogaf.unip.sym$V2 == gogaf.unip.sym$V3),]

head(sort(table(gogaf.unip.sym$V3), decreasing = T), 25)
head(sort(table(gogaf.unip.sym.u$V3), decreasing = T), 25)

hsapuniprot <- fread('../data/hsap.uniprotsymbol.map')

hgncsyms <- fread('../data/public/GOgaf/hgnc.txt', header = T)
setnames(hgncsyms, c('id','apsym','name','prevsym','alias'))
hsap.genes <- readLines('../data/hsap.genes')

toppairs <- fread('../data/toperc/sumlogpftercgt15.bestbitscorepclt10.tsv', header = T)
ids1 <- toppairs$g1 %in% gogaf.bygene$V3
ids2 <- toppairs$g2 %in% gogaf.bygene$V3
ids <- ids1&ids2
View(toppairs[!ids,])

#unmapgenes <- unique(c(toppairs$g1[!ids1], toppairs$g2[!ids2]))
unmapgenes <- hsap.genes[!hsap.genes %in% gogaf.bygene$V3]

map1 <- filter(hgncsyms, prevsym %in% unmapgenes, apsym %in% gogaf.bygene$V3) %>%
     dplyr::select(c(apsym,prevsym)) %>% unique()
map2 <- filter(hgncsyms, apsym %in% unmapgenes, prevsym %in% gogaf.bygene$V3) %>%
     dplyr::select(c(apsym,prevsym)) %>% unique()
mapvec <- c(map1$apsym,map2$prevsym)
names(mapvec) <- c(map1$prevsym,map2$apsym)

saveRDS(mapvec, '../data/hsap.unmapgogaf.genes.rds')

toppairs$g1[!ids1] <- mapvec[toppairs$g1[!ids1]]
toppairs$g2[!ids2] <- mapvec[toppairs$g2[!ids2]]
toppairs <- toppairs %>% na.omit()

mgs1k <- c(mgs,sapply(101:1000, function(x){
     #x <- 1
     print(x)
     mgoSim(filter(gogaf, V3 == toppairs$g1[x])$V5,
            filter(gogaf, V3 == toppairs$g2[x])$V5)
}))

hist(mgs, breaks = 20, probability = T)
hist(mgs1k, breaks = 20, probability = T, col = rgb(1,0,0, alpha = 0.25), add = T)

saveRDS(mgs1k, file = '../data/geneSimBP/top1k.genesim.bp.rds')
     
indf <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T) %>%
     dplyr::select(c('g1','g2','pair','interc')) %>%
     na.omit()

toppairs$pair <- calcpairid(toppairs$g1, toppairs$g2, mat = pairmap)

toppairs <- toppairs %>%
     mutate(interc = qnorm(-1*value, lower.tail = F, log.p = T)) %>%
     dplyr::select(c('g1','g2','pair','interc'))
write.table(toppairs, file = '../data/toperc/top.sumlogpgt15.bestbitscorepclt10.tsv',
            quote = F, sep = '\t', row.names = F)


mapvec <- readRDS('../data/hsap.unmapgogaf.genes')

# d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
# mgoSim(filter(gogaf,V3 == 'CFAP57')$V5,
#        filter(gogaf,V3 == 'DNAH10')$V5)
# 
# head(gogaf.bygene)
# goSim('GO:0002576','GO:0003674')
# 
# 
library(topGO)
library(GOSemSim)
library(org.Hs.eg.db)
# 
# xx <- as.list(org.Hs.egALIAS2EG)
# 
# toppairs <- fread('../data/toperc/sumlogpftercgt15.bestbitscorepclt10.tsv', header = T)
# bgdpairs <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T)
# 
# hsap.genes <- readLines('../data/hsap.genes')
# 
# xx.hsap <- xx[names(xx) %in% hsap.genes]
# 
# geneSim("1","9", ont = 'BP')
# geneSim("1","1982", ont = 'BP')
# geneSim("1","6530", ont = 'BP')
# geneSim("1","10991", ont = 'BP')
# 
# xx.hsap.multi <- xx.hsap[sapply(xx.hsap, length) > 1]
# sort(names(xx.hsap.multi))
