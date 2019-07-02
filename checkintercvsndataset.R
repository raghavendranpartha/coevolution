# meanpfterc <- readRDS('../data/allmeanp.fterc.rds')
# meanpfterc[meanpfterc==0] <- NA
# meanpfterc <- qnorm(meanpfterc,lower.tail = F)
# diag(meanpfterc) <- 1
# saveRDS(meanpfterc, file = '~/Documents/coevolution/data/allmeanpfterc.rds')

#intfterc.ndata <- readRDS('../data/allsum.fterc.rds')
# intfterc.ndata[intfterc.ndata == 0] <- NA
# diag(intfterc.ndata) <- NA
# saveRDS(intfterc.ndata, file = '../data/allsum.fterc.rds')


#badblastpairs <- fread('../data/allgeneshg19.pairbitscore.gteq0p1.tsv')
# intfterc <- readRDS('../data/allsum.fterc.rds')
# ndatafterc <- readRDS('../data/all.ndatasetsfterc.rds')
# dim(ndatafterc)
# ndatafterc[1:5,1:5]
# diag(ndatafterc) <- NA
# sum(ndatafterc == 5, na.rm = T)
# 
# 
# #intfterc.ndata['ABCA10','A2M']
# 
# 
# sum(is.na(intfterc.ndata))
# intmeanfterc.ndata <- intfterc.ndata/ndatafterc

# a <- matrix(c(2,1,2,4), nrow = 2)
# b <- matrix(c(1,4,2,2), nrow = 2)
# a/b

# p1 <- 0.01
# p2 <- 0.05
# curr <- -log(p1)-log(p2)
# curr2 <- -log(mean(c(p1,p2)))
# curr3 <- curr/2


# toppairs <- fread('../data/toperc/top.sumlogpgt15.bestbitscorepclt10.tsv')
# onembgd <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv')
# toppairs$ndata <- factor(calcpairid(toppairs$g1, toppairs$g2, ndatafterc))
# onembgd$ndata <- factor(calcpairid(onembgd$g1, onembgd$g2, ndatafterc))
# table(toppairs$ndata)/sum(table(toppairs$ndata))
# table(onembgd$ndata)/sum(table(onembgd$ndata))
# 
# toppairs <- filter(toppairs, !pair %in% badblastpairs$pair)
# g <- ggplot(toppairs, aes(x=ndata,y=interc))+geom_boxplot()
# print(g)



#intfterc.ndata <- meanpfterc
# sum(rownames(intfterc.ndata) == rownames(ndatafterc))
# ndatafterc[upper.tri(ndatafterc)] <- NA


# nforp <- 2000000
# inds5 <- which(ndatafterc == 5, arr.ind = T)
# inds4 <- which(ndatafterc == 4, arr.ind = T)
# inds3 <- which(ndatafterc == 3, arr.ind = T)
# inds2 <- which(ndatafterc == 2, arr.ind = T)
# inds1 <- which(ndatafterc == 1, arr.ind = T)
# inds3 <- inds3[sample(1:dim(inds3)[1], nforp),]
# inds2 <- inds2[sample(1:dim(inds2)[1], nforp),]
# inds1 <- inds1[sample(1:dim(inds1)[1], nforp),]
#hist(intfterc.ndata[inds5], breaks = 25, col = rgb(1,0,0,0.25), probability = T)
#hist(intfterc.ndata[inds4], breaks = 25, col = rgb(0,0,1,0.25), probability = T, add = T)

bgdbig <- readRDS('../data/toperc/bgdbig.rds')

par(mfrow = c(1,2))
boxplot(
     filter(bgdbig, ndataset == 1)$interc,
     filter(bgdbig, ndataset == 2)$interc,
     filter(bgdbig, ndataset == 3)$interc,
     filter(bgdbig, ndataset == 4)$interc, 
     filter(bgdbig, ndataset == 5)$interc,
     xlab='#datasets',ylab = 'Integrated ERC', main = 'Integrated sum FTerc vs #datasets')
abline(h=median(filter(bgdbig, ndataset == 1)$interc, na.rm = T), lty = 2)
boxplot(
     filter(bgdbig, ndataset == 1)$merc,
     filter(bgdbig, ndataset == 2)$merc,
     filter(bgdbig, ndataset == 3)$merc,
     filter(bgdbig, ndataset == 4)$merc, 
     filter(bgdbig, ndataset == 5)$merc,
     xlab='#datasets',ylab = 'Mammal ERC', main = 'Mammal FTerc vs #datasets')
abline(h=median(filter(bgdbig, ndataset == 1)$merc, na.rm = T), lty = 2)
#abline(h=median(intfterc.ndata[inds4], na.rm = T))

for(ii in c(1:5)) {print(median(merc.sym[get(paste0('inds',ii))], na.rm = T))}
for(ii in c(1:5)) {print(median(intfterc[get(paste0('inds',ii))], na.rm = T))}



wilcox.test(intfterc[inds5], intfterc[inds1], alternative = 'g')
wilcox.test(intfterc[inds4], intfterc[inds1], alternative = 'g')
wilcox.test(intfterc[inds3], intfterc[inds1], alternative = 'g')
wilcox.test(intfterc[inds2], intfterc[inds1], alternative = 'g')


boxplot(
     intmeanfterc.ndata[inds1],
     intmeanfterc.ndata[inds2],
     intmeanfterc.ndata[inds3],
     intmeanfterc.ndata[inds4], 
     intmeanfterc.ndata[inds5],
     xlab='#datasets',ylab = 'Integrated ERC', main = 'Integrated mean FTerc vs #datasets')
abline(h=median(intmeanfterc.ndata[inds1], na.rm = T))
abline(h=median(intmeanfterc.ndata[inds4], na.rm = T))
abline(h=1,lty=2)

#for z-score of sumlogpfterc
par(mfrow = c(1,1))
boxplot(
     qnorm(-1*intfterc.ndata[inds1], lower.tail = F, log.p = T),
     qnorm(-1*intfterc.ndata[inds2], lower.tail = F, log.p = T),
     qnorm(-1*intfterc.ndata[inds3], lower.tail = F, log.p = T),
     qnorm(-1*intfterc.ndata[inds4], lower.tail = F, log.p = T), 
     qnorm(-1*intfterc.ndata[inds5], lower.tail = F, log.p = T),
     xlab='#datasets',ylab = 'Integrated ERC', main = 'Integrated z-score from sum logpFTerc vs #datasets')
abline(h=median(qnorm(-1*intfterc.ndata[inds1], lower.tail = F, log.p = T), na.rm = T))
abline(h=median(qnorm(-1*intfterc.ndata[inds4], lower.tail = F, log.p = T), na.rm = T))
abline(h=2.5,lty=2)
