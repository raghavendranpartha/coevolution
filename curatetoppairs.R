require(dplyr)
require(data.table)
require(reshape2)
require(ggplot2)
require(tidyr)

# blastpair <- fread('../data/allgeneshg19.pairbitscore.tsv', header = T)
# blastpair$pair <- calcpairid(blastpair$g1,blastpair$g2,pairmap)
# badblastpairs <- filter(blastpair, bestbitscorepc >= 0.1)

zfg <- sort(readLines('../data/public/zincfinger.gene'))
org <- sort(readLines('../data/public/olfactoryreceptor.gene'))
orgzfg <- sort(c(org,zfg))
# org.c <- paste(t(combn(org,2))[,1], t(combn(org,2))[,2], sep = ',')
# zfg.c <- paste(t(combn(zfg,2))[,1], t(combn(zfg,2))[,2], sep = ',')
#orgzfg.c <- paste(t(combn(orgzfg,2))[,1], t(combn(orgzfg,2))[,2], sep = ',')

top.df <- melt(intfterc) %>% filter(value >= 5) %>%
     mutate(Var1 = as.character(Var1)) %>% mutate(Var2 = as.character(Var2)) %>%
     filter(!Var1 %in% orgzfg, !Var2 %in% orgzfg) %>%
     mutate(pair = calcpairid(Var1,Var2,pairmap)) %>% group_by(pair) %>%
     summarise(Var1 = Var1[1], Var2 = Var2[1], value = value[1]) %>%
     arrange(desc(value))
setnames(top.df, c('pair','g1','g2','interc'))
write.table(top.df, file = '../data/toperc/top.sumfterc5.bestbitscorepclt10.tsv',
            quote = F, sep = '\t', row.names = F)

top.df <- fread('../data/toperc/top.sumfterc5.bestbitscorepclt10.tsv')
bgdbig <- readRDS('../data/toperc/bgdbig.rds')
bgdbig <- bgdbig[!is.na(bgdbig$interc),]
top.df$ndataset <- calcpairid(top.df$g1, top.df$g2, ndatafterc)
ntop <- 100000
control.top.df <- nds.matched.controlpairs(top.df[1:ntop,]$ndataset, nPerLevelmult = 5)
sum(is.na(control.top.df$interc))
write.table(control.top.df, file = '../data/toperc/controlndatasetsmatched.top.sumfterc5.bestbitscorepclt10.tsv',
            quote = F, sep = '\t', row.names = F)


top.df <- melt(intfterc) %>% filter(value <= -5) %>%
     mutate(Var1 = as.character(Var1)) %>% mutate(Var2 = as.character(Var2)) %>%
     filter(!Var1 %in% orgzfg, !Var2 %in% orgzfg) %>%
     mutate(pair = calcpairid(Var1,Var2,pairmap)) %>% group_by(pair) %>%
     summarise(Var1 = Var1[1], Var2 = Var2[1], value = value[1]) %>%
     arrange(desc(value))
setnames(top.df, c('pair','g1','g2','interc'))
write.table(top.df, file = '../data/toperc/negativetop.sumftercm5.bestbitscorepclt10.tsv',
            quote = F, sep = '\t', row.names = F)
     
#ignore the rest


# alllogpfterc[upper.tri(alllogpfterc)] <- NA
# diag(alllogpfterc) <- NA
# alllogpfterc.df <- melt(alllogpfterc) %>% filter(value >= 5) %>%
#      mutate(Var1 = as.character(Var1)) %>% mutate(Var2 = as.character(Var2)) %>%
#      mutate(pair = calcpairid(Var1,Var2)) 
# rm('alllogpfterc')
# alllogpfterc.df <- left_join(alllogpfterc.df, select(blastpair, c('pair','bestbitscorepc'))) %>% arrange(desc(value))
# alllogpfterc.df$bestbitscorepc[is.na(alllogpfterc.df$bestbitscorepc)] <- 0.000001
# topscoringpairs <- filter(alllogpfterc.df,!Var1 %in% orgzfg,!Var2 %in% orgzfg, bestbitscorepc < 0.1)
# write.table(topscoringpairs, file = '../data/toperc/sumlogpftercgt5.bestbitscorepclt10.tsv',
#             quote = F, sep = '\t', row.names = F)

g <- ggplot(topscoringpairs, aes(x=value, y=bestbitscorepc))+
     geom_point()+geom_smooth(method='lm')
print(g)

with(filter(alllogpfterc.df,!g1 %in% orgzfg,!g2 %in% orgzfg, bestbitscorepc < 0.1),
     head(sort(table(c(g1,g2)), decreasing = T), 100))






npairs <- 1200000
set.seed(2)
allpairs <- cbind(sample(hsap.genes,npairs,replace = T),sample(hsap.genes,npairs,replace = T))
allpairs <- allpairs[(allpairs[,1] != allpairs[,2]),]
res <- paste(allpairs[,1], allpairs[,2], sep = ',')
resrev <- paste(allpairs[,2], allpairs[,1], sep = ',')
badpairs <- unique(c(intersect(blastpair$pair, res),
              intersect(blastpair$pair, resrev)))
badinds <- match(badpairs, res)
badindsrev <- match(badpairs, resrev)
allpairs <- allpairs[-c(badinds[!is.na(badinds)],badindsrev[!is.na(badindsrev)]),]
resfinal <- unique(res[-c(badinds[!is.na(badinds)],badindsrev[!is.na(badindsrev)])])
rm(list=c('res','resrev','badpairs','badinds','badindsrev','allpairs'))
resfinal.mat <- cbind(sapply(strsplit(resfinal,','),'[[',1),
                      sapply(strsplit(resfinal,','),'[[',2))
rm('resfinal')
rowinds <- match(resfinal.mat[,1], rownames(alllogpfterc))
colinds <- match(resfinal.mat[,2], colnames(alllogpfterc))
res.sumlogpfterc <- alllogpfterc[rowinds + nrow(alllogpfterc) * (colinds - 1)]


nrw <- nrow(blastpair)
print(Sys.time())
rowIndices <- match(blastpair$g1[1:nrw], rownames(alllogpfterc))
colIndices <- match(blastpair$g2[1:nrw], colnames(alllogpfterc))
blastpair$sum.logpfterc <- alllogpfterc[rowIndices + nrow(alllogpfterc) * (colIndices - 1)]
print(Sys.time())

boxplot(res.sumlogpfterc, log = 'y')

g <- ggplot(blastpair, aes(x = bestbitscorepc, sum.logpfterc))+geom_point()+geom_smooth(method = 'lm')+coord_trans(y='log10')+theme_bw()
print(g)


resfinal.df <- data.frame(sum.logpfterc = res.sumlogpfterc,
                          bestbitscorepc = 0.0001, stringsAsFactors = F)
plotdf <- rbind(select(blastpair,c('sum.logpfterc','bestbitscorepc')),
                resfinal.df)

nbreaks <- 19
qq=quantile(plotdf$bestbitscorepc,seq(0,nbreaks,1)/nbreaks)
qqdiff=diff(qq)
breaks=qq[1:nbreaks]+qqdiff/2
rr=quantile(plotdf$bestbitscorepc, c(0.0001, 0.99))
breaks=round(breaks,5)
breaks <- seq(range(blastpair$bestbitscorepc)[1],range(blastpair$bestbitscorepc)[2], length.out = nbreaks)

breaks=c(0,round(breaks,5))
breaks <- c(0,0.0001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.4,0.8,1)
nbreaks <- length(breaks)-1

par(mar = c(10,5,3,1))
cutres<-cut(plotdf$bestbitscorepc,breaks = breaks)
cutres_tt <- table(cutres)
boxplot((plotdf$sum.logpfterc)~ cutres, xlab = "", ylab = "sum.logp.fterc", outline=F,  log="y", las=2, srt = 45, xaxt='n')
abline(h=median(resfinal.df$sum.logpfterc, na.rm = T), lty = 2, col = 'red')
axis(1, at=1:(nbreaks), labels = F)
# text(1:(nbreaks-1), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), srt = 45, adj = 1,
#      labels = levels(cutres), xpd = TRUE, cex = 1.5)
text(1:(nbreaks), y =1e-09, srt = 45, adj = 1,
     labels = levels(cutres), xpd = TRUE, cex = 1.25)
text(1:(nbreaks), y =1e-07, srt = 30, adj = 0.5,
     labels = cutres_tt, xpd = TRUE, cex = 1.0)
title("Integrated ERC Score vs Relative Blast bit score")



#zinc fingers
amingozfg <- fread('../data/public/amigozincfinger.txt', header = F) %>%
     filter(grepl('zinc finger protein', V3, ignore.case = T))
zfg1 <- intersect(unique(amingozfg$V2),hsap.genes)
zfg2 <- grep('^ZNF',hsap.genes, value = T)
zfg3 <- grep('^ZFP',hsap.genes, value = T)

zfg <- unique(c(zfg1,zfg2,zfg3))

orgenes <- sort(grep('^OR',hsap.genes,value = T))
amigoor <- fread('../data/public/amigoolfactoryreceptor.txt', header = F) %>%
     filter(grepl('olfactory receptor', V2, ignore.case = T))
org1 <- intersect(unique(amigoor$V3),hsap.genes)
setdiff(orgenes,org1)
annots <- fread('~/Documents/createAlns/data/genedescriptions/hgTables', header = T) %>%
     select(c(7,8)) %>% setnames(c('gene','summary')) %>% unique()
org2 <- filter(annots, gene %in% setdiff(orgenes, org1),grepl('olfactory receptor',summary))$gene
org <- unique(c(org1,org2))

write(zfg, file = '../data/public/zincfinger.gene')
write(org, file = '../data/public/olfactoryreceptor.gene')

#for zinc finger proteins     
zfg.blastpair <- filter(blastpair, g1 %in% zfg, g2 %in% zfg)
#for olfactory receptor proteins
zfg.blastpair <- filter(blastpair, g1 %in% org, g2 %in% org)

plotdf <- rbind(select(zfg.blastpair,c('sum.logpfterc','bestbitscorepc')),
                resfinal.df)
par(mar = c(10,5,3,1))
breaks <- c(0,0.0001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.4,0.8,1)
nbreaks <- length(breaks)-1
cutres<-cut(plotdf$bestbitscorepc,breaks = breaks)
cutres_tt <- table(cutres)
boxplot((plotdf$sum.logpfterc)~ cutres, xlab = "", ylab = "sum.logp.fterc", outline=F,  log="y", las=2, srt = 45, xaxt='n')
abline(h=median(resfinal.df$sum.logpfterc, na.rm = T), lty = 2, col = 'red')
axis(1, at=1:(nbreaks), labels = F)
# text(1:(nbreaks-1), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), srt = 45, adj = 1,
#      labels = levels(cutres), xpd = TRUE, cex = 1.5)
text(1:(nbreaks), y =1e-09, srt = 45, adj = 1,
     labels = levels(cutres), xpd = TRUE, cex = 1.25)
text(1:(nbreaks), y =1e-07, srt = 30, adj = 0.5,
     labels = cutres_tt, xpd = TRUE, cex = 1.0)
#title("Zinc Finger Proteins\nIntegrated ERC Score vs Relative Blast bit score")
title("Olfactory receptor Proteins\nIntegrated ERC Score vs Relative Blast bit score")


zfg <- sort(readLines('../data/public/zincfinger.gene'))
org <- sort(readLines('../data/public/olfactoryreceptor.gene'))

org.c <- paste(t(combn(org,2))[,1], t(combn(org,2))[,2], sep = ',')
zfg.c <- paste(t(combn(zfg,2))[,1], t(combn(zfg,2))[,2], sep = ',')


zfg.pair <- expand.grid(zfg,zfg) %>%mutate(Var1 = as.character(Var1)) %>% mutate(Var2 = as.character(Var2)) %>%
     mutate(pair = paste(sort(c(Var1,Var2)),collapse=',')) %>% unique()
