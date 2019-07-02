require(dplyr)
require(data.table)

# dgg.np <- readLines('../data/public/omim/omim.npreid.2015.list') 
# dgg.names <- readLines('../data/public/omim/omim.npreid.2015.disease')
# dgg.np.g <- dgg.np[seq(2,1290,by=4)]
# dgg.np.g.use <- strsplit(dgg.np.g,",",fixed = F)
# dgg.u <- cleanedcomplex(dgg.np.g.use, dgg.names)

#dgg.np.g.use[dgg.names == 'Monilethrix']
#intersect(dgg.np.g.use[dgg.names == 'Monilethrix'],hsap.genes.useforsim)

#saveRDS(dgg.u, file = '../data/public/omimcplxs.cleaned.rds')
dgg.u <- readRDS('../data/public/omimcplxs.cleaned.rds')

#dgg.u$df$ndata <- calcpairid(dgg.u$df$V1,dgg.u$df$V2,ndatafterc)
View(dgg.u$df)

dgg.df.summ <- group_by(dgg.u$df, cplxid) %>%
     summarise(meanintscore = mean(interc, na.rm = T),
               meanmerc = mean(merc, na.rm = T)) 

randfiles <- dir('../data/randoms/omimcplxs.cleaned/')
randcplxdfsumm <- do.call('rbind',lapply(randfiles, function(x){
     fread(file.path('../data/randoms/omimcplxs.cleaned',x))
}))
#head(randcplxdfsumm)
# randcplxdfsumm.mean <- group_by(randcplxdfsumm, cplxid) %>%
#      summarise(meanofmeans = mean(medianscore, na.rm = T))

dgg.df.summ$pvint <- sapply(1:nrow(dgg.df.summ), function(x){
     mean(filter(dgg.df.summ, cplxid == x)$meanintscore<=filter(randcplxdfsumm,cplxid == x)$meanintscore,
          na.rm = T)
})
dgg.df.summ$pvmamm <- sapply(1:nrow(dgg.df.summ), function(x){
     mean(filter(dgg.df.summ, cplxid == x)$meanmerc<=filter(randcplxdfsumm,cplxid == x)$meanmerc,
          na.rm = T)
})
dgg.df.summ$pvint[dgg.df.summ$pvint == 0] <- 1e-06
dgg.df.summ$pvmamm[dgg.df.summ$pvmamm == 0] <- 1e-06
dgg.df.summ$names <- names(dgg.u$cplx)
dgg.df.summ$ngenes <- sapply(dgg.u$cplx, length)
dgg.df.summ$qvint <- p.adjust(dgg.df.summ$pvint, method = 'BH')
dgg.df.summ$qvmamm <- p.adjust(dgg.df.summ$pvmamm, method = 'BH')

head(dgg.df.summ)
dgg.df.summ <- select(dgg.df.summ,c(1,6,7,2:5,8,9))
setnames(dgg.df.summ, c('Complex.ID','Disease group','nGenes','Mean.Integrated.ERC','Mean.Mammal.ERC','P.value.Integrated','P.value.Mammal','Q.value.Integrated',
                         'Q.value.Mammal'))

write.table(arrange(dgg.df.summ,Q.value.Integrated),'../data/tables/omim.all.tsv', 
            quote = F, sep = '\t', row.names = F)

dgg.df.summ <- fread('../data/tables/omim.all.tsv')

sum(dgg.df.summ$Q.value.Integrated < 0.05, na.rm = T)
sum(dgg.df.summ$Q.value.Mammal < 0.05, na.rm = T)

dgg.df.summ.forroc <- select(dgg.df.summ, c(4,5))
setnames(dgg.df.summ.forroc, c('interc','merc'))
setnames(randcplxdfsumm,c('cplxid','dataset','interc','merc'))
# onembgd <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T) %>%
#      filter(!pair %in% dgg.u$df$pair)
head(dgg.u$df)

#pdf('../figures/aurocs/omimvsrandom.pdf',width = 5,height = 5,family = 'FreeSans')
#dev.off()
#runroc.intvsmamm(dgg.u$df, onembgd,plot = T, titl = 'OMIM DGG PPIs vs Random PPIs')

#calc full auc
#runroc.intvsmamm(dgg.df.summ.forroc, filter(randcplxdfsumm, dataset <= 100000),plot = T, titl = '')

pdf('../figures/pvals/omimvsrandomPvalsauroc.pdf',width = 10,height = 5,family = 'FreeSans')
par(mfrow=c(1,2))
hist(dgg.df.summ$P.value.Integrated, breaks = 25, probability = T, xlab = 'p-value', col = rgb(1,0,0,0.3), main = '',
     cex.lab = 1.25)
hist(dgg.df.summ$P.value.Mammal, breaks = 25, probability = T, col = rgb(0,0,0,0.3), add = T)
title('A',adj = 0, font.main = 1,cex.main = 1.5)
legend('right',c('Integrated','Mammal'),title = 'ERC',col = c('red','black'),lty = 1,bty='n',cex = 1.25)
#runroc.intvsmamm(dgg.df.summ.forroc, filter(randcplxdfsumm, dataset <= 10),plot = T, titl = 'B')
runpr.intvsmamm(dgg.df.summ.forroc, randcplxdfsumm,plot = T, titl = 'B',ypts = c(1,5,10,50,250,1000))
dev.off()

#IGNOre

qqplot(x = randcplxdfsumm.mean$meanofmeans, y = dgg.df.summ$medianscore, xlab = 'Mean of means in 100k random complexes', ylab = 'OMIM DGGs mean score',
       xlim = round(range(dgg.df.summ$medianscore, na.rm = T),3))
abline(a=0,b=1,lty=2)

nrandoms <- 500
randcplxs <- do.call('rbind',lapply(1:nrandoms, function(y){
     print(y)
     createRandomCplx.df(hsap.genes.useforsim, dgg) %>% mutate(dataset = y)
}))
randcplxs$pair <- with(randcplxs,calcpairid(V1,V2,pairmap))
randcplxs$sumlogpfterc <- with(randcplxs,calcpairid(V1,V2,alllogpfterc))
randcplxs$merc <- with(randcplxs,calcpairid(V1,V2,merc.sym))

dgg.df$lbl <- 1
randcplxs$lbl <- 0

dgg.df$interc <- qnorm(-1*dgg.df$sumlogpfterc, lower.tail = F, log.p = T)
randcplxs$interc <- qnorm(-1*randcplxs$sumlogpfterc, lower.tail = F, log.p = T)

nball <- runNaiveBayes.intvsmamm(dgg.df, randcplxs,plot=T,titl = 'OMIM DGGs')


randcplxs.summ <- group_by(randcplxs, cplxid,dataset) %>%
     summarise(medianscore = mean(interc, na.rm = T))
dgg.df.summ$complex <- 'OMIM DGG'
randcplxdfsumm$complex <- 'Random'
plotdf <- bind_rows(dgg.df.summ, dplyr::select(randcplxdfsumm,c(1,3,4)))
plotdf$complex <- factor(plotdf$complex, levels = c('Random','OMIM DGG'))




sum(dgg.df.summ$qv < 0.05, na.rm = T)

g <- ggplot(plotdf, aes(x = complex, y=medianscore, fill = complex))+
     geom_violin()+geom_boxplot(width=0.2)+theme_bw()+
     xlab('')+ylab('Mean Integrated Score')+ggtitle('OMIM DGGs')+
     theme(axis.text=element_text(size = 18, face='bold'),
           axis.title=element_text(size = 20, face='bold'),
           plot.title=element_text(size = 22, face='bold'),
           strip.text = element_text(size = 18, face = 'bold'))
print(g)


