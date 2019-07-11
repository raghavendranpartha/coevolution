library(RERconverge)
require(dplyr)
require(data.table)
require(ggplot2)

mtre <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
vtre <- readRDS('~/Documents/erc/data/saves/vert39.trees.rds')
dtre <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')
wtre <- readRDS('~/Documents/erc/data/saves/worms17.trees.rds')
ytre <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.rds')

source('../code/funcsERCfromProjections.R')

hsap.genes <- readLines('../data/hsap.genes')

intfterc <- readRDS('../data/allsum.fterc.rds')
merc.sym <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.symm.rds')
pairmap <- readRDS('../data/pairmapmat.rds')
zfg <- sort(readLines('../data/public/zincfinger.gene'))
org <- sort(readLines('../data/public/olfactoryreceptor.gene'))
orgzfg <- sort(c(org,zfg))
hsap.genes.useforsim <- setdiff(hsap.genes,orgzfg)
badblastpairs <- fread('../data/allgeneshg19.pairbitscore.gteq0p1.tsv', header = T)

ndatafterc <- readRDS('../data/all.ndatasetsfterc.rds')
#alllogpfterc <- readRDS('~/Documents/coevolution/data/allsum.logpfterc.rds')




blastpair <- fread('../data/allgeneshg19.pairbitscore.tsv', header = T)
blastpair$pair <- calcpairid(blastpair$g1,blastpair$g2,pairmap)
badblastpairs <- filter(blastpair, bestbitscorepc >= 0.1)
write.table(badblastpairs, file = '../data/allgeneshg19.pairbitscore.gteq0p1.tsv', 
            quote = F, sep = '\t', row.names = F)
rm('blastpair')

merc <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.rds')
verc <- readRDS('~/Documents/erc/data/saves/vert39.trees.fterc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.fterc.rds')
werc <- readRDS('~/Documents/erc/data/saves/worms17.trees.fterc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.fterc.rds')
mapvecs <- readRDS('../data/flywormyeastmap.vec.rds')

dercmap <- mappedERC(derc, mapvecs[[1]])
wercmap <- mappedERC(werc, mapvecs[[2]])
yercmap <- mappedERC(yerc, mapvecs[[3]])

omerc <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.erc.rds')


par(mar = c(1,1,0,0))
plot(ytre$masterTree, font = 2,cex = 0.75)

mrers <- readRDS('~/Documents/rermethods/data/mamm63nt.trees.scaledrers.sqrt.wt.rds')

wormogwbgn <- fread('../data/cele.OgWbgnTxID.map')

flymap <- fread('../data/dmel.hsap.geneIDsUniProtIDs.map')
yeastmap <- fread('../data/hsap.scer.geneIDsUniProtIDs.map')
wormmap <- fread('../data/cele.OgwbgnUniprotSymbol.map')
wormmap.vec <- wormmap$og
names(wormmap.vec) <- wormmap$symbol
flymap.vec <- flymap$og
names(flymap.vec) <- flymap$hsapsymbol
yeastmap.vec <- yeastmap$scersymbol
names(yeastmap.vec) <- yeastmap$hsapsymbol

werc[is.infinite(werc)] <- NA
celeids <- which(rownames(werc) %in% wormogwbgn$og)

wbgnmap <- wormogwbgn$V2
names(wbgnmap) <- wormogwbgn$og

rownames(werc)[celeids] <- wbgnmap[rownames(werc)[celeids]]
colnames(werc)[celeids] <- wbgnmap[colnames(werc)[celeids]]

saveRDS(werc, file = '../data/erc/worms17.erc.wbgnids.rds')

View(as.data.frame(rownames(werc)))

#nvals <- 10000
#hist(mtre$paths, breaks = 25, probability = T)
qqm = quantile(mtre$paths, probs = seq(0,1,length.out = 101), na.rm = T)
qqv = quantile(vtre$paths, probs = seq(0,1,length.out = 101), na.rm = T)
qqd = quantile(dtre$paths, probs = seq(0,1,length.out = 101), na.rm = T)
qqw = quantile(wtre$paths, probs = seq(0,1,length.out = 101), na.rm = T)
qqy = quantile(ytre$paths, probs = seq(0,1,length.out = 101), na.rm = T)

lins <- c('vert','dmel','worm','scer')
cols <- c('red','blue','magenta','green')
plot(seq(0,1,length.out = 101), qqm, log = 'y',type = 'l', lwd = 2, ylim = c(1e-06,100), xlab = 'probability',
     ylab = 'quantiles')
lines(seq(0,1,length.out = 101), qqv, lwd = 2, col = 'red')
lines(seq(0,1,length.out = 101), qqd, lwd = 2, col = 'blue')
lines(seq(0,1,length.out = 101), qqw, lwd = 2, col = 'magenta')
lines(seq(0,1,length.out = 101), qqy, lwd = 2, col = 'green')
legend('topleft',c('mamm',lins),col = c('black',cols),lty = 1,lwd=2)