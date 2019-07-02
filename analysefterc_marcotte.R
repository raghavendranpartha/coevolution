#marc.cplx <- readRDS('../data/marcotte981complexes.rds')

require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)
require(beanplot)
library(e1071)
require(PRROC)
require(precrec)

dat <- fread('../data/Marcotte_981complexes.tsv', sep='\t') 
cplx <- strsplit(dat$GeneName,';')
cplx.u <- cleanedcomplex(cplx)
saveRDS(cplx.u, file = '../data/public/marcottecplxs.cleaned.rds')

genepresenceint <- readRDS('../data/genepresenceasinteger.rds')
marc.u <- readRDS('../data/public/marcottecplxs.cleaned.rds')

#randcplxmarc.u <- createRandomCplx.df(hsap.genes.useforsim, marc.u$cplx)

marc.df.summ <- group_by(marc.u$df, cplxid) %>%
     summarise(meanintscore = mean(interc, na.rm = T),
               meanmerc = mean(merc, na.rm = T)) 

randfiles <- dir('../data/randoms/marcottecplxs.cleaned/')
randcplxdfsumm <- do.call('rbind',lapply(randfiles, function(x){
     fread(file.path('../data/randoms/marcottecplxs.cleaned',x))
}))
#head(randcplxdfsumm)
# randcplxdfsumm.mean <- group_by(randcplxdfsumm, cplxid) %>%
#      summarise(meanofmeans = mean(medianscore, na.rm = T))

marc.df.summ$pvint <- sapply(1:nrow(marc.df.summ), function(x){
     mean(filter(marc.df.summ, cplxid == x)$meanintscore<=filter(randcplxdfsumm,cplxid == x)$meanintscore,
          na.rm = T)
})
marc.df.summ$pvmamm <- sapply(1:nrow(marc.df.summ), function(x){
     mean(filter(marc.df.summ, cplxid == x)$meanmerc<=filter(randcplxdfsumm,cplxid == x)$meanmerc,
          na.rm = T)
})
marc.df.summ$pvint[marc.df.summ$pvint == 0] <- 1e-06
marc.df.summ$pvmamm[marc.df.summ$pvmamm == 0] <- 1e-06
#marc.df.summ$names <- names(marc.u$cplx)
marc.df.summ$ngenes <- sapply(marc.u$cplx, length)
marc.df.summ$qvint <- p.adjust(marc.df.summ$pvint, method = 'BH')
marc.df.summ$qvmamm <- p.adjust(marc.df.summ$pvmamm, method = 'BH')

marc.df.summ <- select(marc.df.summ, c(1,6,2:5,7,8))
head(marc.df.summ)
setnames(marc.df.summ, c('Complex.ID','nGenes','Mean.Integrated.ERC','Mean.Mammal.ERC','P.value.Integrated','P.value.Mammal','Q.value.Integrated',
                         'Q.value.Mammal'))

write.table(arrange(marc.df.summ,Q.value.Integrated),'../data/tables/marcotte.all.tsv', 
            quote = F, sep = '\t', row.names = F)

marc.df.summ <- fread('../data/tables/marcotte.all.tsv')

marc.df.summ.forroc <- select(marc.df.summ, c(3,4))
setnames(marc.df.summ.forroc, c('interc','merc'))
setnames(randcplxdfsumm,c('cplxid','dataset','interc','merc'))
head(marc.u$df)
#just to calculate the full aucs
runroc.intvsmamm(marc.df.summ.forroc, filter(randcplxdfsumm, dataset <= 100000),plot = T, titl = 'Metazoan complexes vs Random complexes')


#runroc.intvsmamm(marc.u$df, onembgd,plot = T, titl = 'Marcotte complex PPIs vs Random PPIs')

pdf('../figures/pvals/marcottevsrandomPvalsAUROC.pdf',width =10,height = 5,family = 'FreeSans')
par(mfrow=c(1,2))
hist(marc.df.summ$P.value.Integrated, breaks = 25, probability = T, xlab = 'p-value', col = rgb(1,0,0,0.3), main = '',cex.lab = 1.25)
hist(marc.df.summ$P.value.Mammal, breaks = 25, probability = T, col = rgb(0,0,0,0.3), add = T)
title('A',font.main = 1, adj = 0, cex.main = 1.5)
legend('right',c('Integrated','Mammal'),title = 'ERC',col = c('red','black'),lty = 1,bty='n',cex = 1.25)
#runroc.intvsmamm(marc.df.summ.forroc, filter(randcplxdfsumm, dataset <= 100),plot = T, titl = 'B')
runpr.intvsmamm(marc.df.summ.forroc, randcplxdfsumm,plot = T, titl = 'B',ypts = c(1,5,10,100,500),ymax=750)
dev.off()




sum(marc.df.summ$Q.value.Integrated < 0.05, na.rm = T)
sum(marc.df.summ$Q.value.Mammal < 0.05, na.rm = T)







#IGNORe
nrandoms <- 10
randcplxs <- do.call('rbind',lapply(1:nrandoms, function(y){
     print(y)
     createRandomCplx.df(hsap.genes.useforsim, cplx.u$cplx) %>% mutate(dataset = y)
}))
randcplxs$pair <- with(randcplxs,calcpairid(V1,V2,pairmap))
badpairids <- randcplxs$pair %in% badblastpairs$pair
randcplxs$interc[badpairids] <- NA
randcplxs$merc[badpairids] <- NA
randcplxs$sumlogpfterc <- with(randcplxs,calcpairid(V1,V2,alllogpfterc))
randcplxs$interc <- qnorm(-1*randcplxs$sumlogpfterc, lower.tail = F, log.p = T)
randcplxs$merc <- with(randcplxs,calcpairid(V1,V2,merc.sym))

cplx.u$cplx.use.pair.df$lbl <- 1
randcplxs$lbl <- 0
nball <- runNaiveBayes.intvsmamm(cplx.u$cplx.use.pair.df, randcplxs,plot=T,titl = 'Marcotte Complex PPIs')

cplx.u.summ <- group_by(cplx.u$df, cplxid) %>%
     summarise(medianscore = mean(interc,na.rm = T))
randcplxs.summ <- group_by(randcplxs, cplxid,dataset) %>%
     summarise(medianscore = mean(interc, na.rm = T))
cplx.u$cplx.use.summ$complex <- 'Marcotte'
randcplxs.summ$complex <- 'Random'
plotdf <- bind_rows(cplx.u$cplx.use.summ, dplyr::select(randcplxs.summ,c(1,3,4)))
plotdf$complex <- factor(plotdf$complex, levels = c('Random','Marcotte'))

g <- ggplot(plotdf, aes(x = complex, y=medianscore, fill = complex))+
     geom_violin()+geom_boxplot(width=0.2)+theme_bw()+
     xlab('')+ylab('Mean Integrated Score')+ggtitle('Marcotte complex median score')+
     theme(axis.text=element_text(size = 18, face='bold'),
           axis.title=element_text(size = 20, face='bold'),
           plot.title=element_text(size = 22, face='bold'),
           strip.text = element_text(size = 18, face = 'bold'))
print(g)
     
marcpvs <- sapply(1:nrow(cplx.hsap.use.df.summ), function(x){
     mean(filter(cplx.hsap.use.df.summ, cplxid == x)$medianscore<=filter(randcplxs.summ,cplxid == x)$medianscore,
          na.rm = T)
})

quantile(marcpvs, probs = seq(0,1,length.out = 10), na.rm = T)

geneomaage <- fread('../data/Marcotte_981complexes.genes.omaage.tsv', header = T)
head(geneomaage)
geneomaagevec <- geneomaage$distOMA
names(geneomaagevec) <- geneomaage$symbol

cplx.hsap.use.age <- sapply(1:length(cplx.hsap.use), function(x){
     mean(geneomaagevec[cplx.hsap.use[[x]]]=='new',na.rm = T)
})
cplx.hsap.use.age <- sapply(1:length(cplx.hsap.use), function(x){
     mean(geneomaagevec[cplx.hsap.use[[x]]],na.rm = T)
})

plotdf2 <- data.frame(pvalue = marcpvs,
                      age = cplx.hsap.use.age, stringsAsFactors = F) %>%
     na.omit()
plotdf2$pvaluelv <- cut(plotdf2$pvalue, breaks = c(0,0.06,1.1), right = F)

g <- ggplot(plotdf2, aes(x=pvaluelv, y=age, fill = pvaluelv))+
     geom_violin()+geom_boxplot(width=0.2)+ggtitle('Marcotte power vs age')+theme_bw()+
     theme(axis.text=element_text(size = 18, face='bold'),
            axis.title=element_text(size = 20, face='bold'),
            plot.title=element_text(size = 22, face='bold'),
            strip.text = element_text(size = 18, face = 'bold'))
print(g)

boxplot(filter(plotdf2,age <= 0.45)$pvalue,plotdf2$pvalue)

# 
# 
# hist(randcplxs.logpfterc.summ$medianscore, breaks = 75, col = rgb(1,0,0,alpha = 0.25), probability = T)
# hist(cplx.hsap.use.logpfterc.summ$medianscore, breaks = 75, col = rgb(0,0,0,alpha = 0.25), add = T, probability = T)
# 
# #individual predictors into NB
# marc.erc <- lapply(1:length(cplx.hsap.use), function(x){
#      print(x)
#      getERCasdf(cplx.hsap.use[[x]]) %>% mutate(cplxid = x)
# })
# marc.erc.df <- spread(do.call('rbind',marc.erc), value = 'value', key = 'dataset')
# marc.erc.df$lbl = 1
# 
# nrandoms <- 100
# #nrandoms <- 5
# randcplxs <- lapply(1:nrandoms, function(y){
#      rand.cplx <- createRandomCplx(hsap.genes, cplx.hsap.use)
#      rand.cplx.erc <- lapply(1:length(rand.cplx), function(x){
#           getERCasdf(rand.cplx[[x]]) %>% mutate(cplxid = x)
#      })
#      rand.cplx.erc.df <- spread(do.call('rbind',rand.cplx.erc), value = 'value', key = 'dataset') %>% dplyr::select(c(2:7))
#      rand.cplx.erc.df$lbl = 0
#      rand.cplx.erc.df
# })
# 
# 
# nball <- lapply(1:nrandoms, function(y){
#      print(y)
#      runNaiveBayes(marc.erc.df, randcplxs[[y]],plot=F)
# })
# 
# nbmarcmmdat <- mmdata(scores = join_scores(lapply(nball, '[[', 1), lapply(nball, '[[', 2)), 
#                       labels = join_labels(rep(lapply(nball, '[[', 3),2)), 
#                       modnames = rep(c('Integrated','Mammal only'), each = nrandoms),
#                       dsids = rep(c(1:nrandoms),2))
# mmcurves <- evalmod(nbmarcmmdat,cb_alpha = 0.05)
# autoplot(mmcurves, "ROC", show_cb = TRUE )+coord_cartesian(xlim = c(0,1))+
#      labs(fill="Method")+
#      theme_bw()+ggtitle('Naive Bayes Prediction: Marcotte Complexes')+
#      theme(axis.text=element_text(size = 18, face='bold'),
#            axis.title=element_text(size = 20, face='bold'),
#            plot.title=element_text(size = 22, face='bold'),
#            strip.text = element_text(size = 18, face = 'bold'))
# 
# autoplot(mmcurves, "PRC", show_cb = TRUE )+coord_cartesian(xlim = c(0,1))+
#      labs(fill="Method")+
#      theme_bw()+ggtitle('Naive Bayes Prediction: Marcotte Complexes')+
#      theme(axis.text=element_text(size = 18, face='bold'),
#            axis.title=element_text(size = 20, face='bold'),
#            plot.title=element_text(size = 22, face='bold'),
#            strip.text = element_text(size = 18, face = 'bold'))
