require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)
source('~/Documents/erc/code/funcsERCfromProjections.R')


allcorumraw <- fread("../data/public/allComplexes.txt", header = T)

corumraw <- fread('../data/public/coreComplexes.txt', header = T)
corum <- corumraw %>% 
     filter(Organism == 'Human') %>%
     dplyr::select(c(1,2,6,13,18))
corumcplx <- strsplit(corum$`subunits(Gene name)`,';')
corumcplx.u <- cleanedcomplex(corumcplx,corum$ComplexName)

saveRDS(corumcplx.u, file = '../data/public/corumcplxs.cleaned.rds')

corum.u <- readRDS(file = '../data/public/corumcplxs.cleaned.rds' )
corum.df.summ <- group_by(corum.u$df, cplxid) %>%
     summarise(meanintscore = mean(interc, na.rm = T),
               meanmerc = mean(merc, na.rm = T)) 

randfiles <- dir('../data/randoms/corumcplxs.cleaned/')
randcplxdfsumm <- do.call('rbind',lapply(randfiles, function(x){
     fread(file.path('../data/randoms/corumcplxs.cleaned',x))
}))
#head(randcplxdfsumm)
# randcplxdfsumm.mean <- group_by(randcplxdfsumm, cplxid) %>%
#      summarise(meanofmeans = mean(medianscore, na.rm = T))

corum.df.summ$pvint <- sapply(1:nrow(corum.df.summ), function(x){
     mean(filter(corum.df.summ, cplxid == x)$meanintscore<=filter(randcplxdfsumm,cplxid == x)$meanintscore,
          na.rm = T)
})
corum.df.summ$pvmamm <- sapply(1:nrow(corum.df.summ), function(x){
     mean(filter(corum.df.summ, cplxid == x)$meanmerc<=filter(randcplxdfsumm,cplxid == x)$meanmerc,
          na.rm = T)
})
corum.df.summ$pvint[corum.df.summ$pvint == 0] <- 1e-06
corum.df.summ$pvmamm[corum.df.summ$pvmamm == 0] <- 1e-06
corum.df.summ$names <- names(corum.u$cplx)
corum.df.summ$ngenes <- sapply(corum.u$cplx, length)
corum.df.summ$qvint <- p.adjust(corum.df.summ$pvint, method = 'BH')
corum.df.summ$qvmamm <- p.adjust(corum.df.summ$pvmamm, method = 'BH')

corum.df.summ <- arrange(select(corum.df.summ, c(1,6,7,2,3,4,5,8,9)),qvint)
setnames(corum.df.summ,c('Complex ID','Name','nGenes','Mean.Integrated.ERC', 'Mean.Mammal.ERC','P.value.Integrated','P.value.Mammal',
                                         'Q.value.Integrated','Q.value.Mammal'))

write.table(corum.df.summ,'../data/tables/corum.all.tsv', 
            quote = F, sep = '\t', row.names = F)

corum.df.summ <- fread('../data/tables/corum.all.tsv', header = T)
corum.df.summ.unq <- group_by(corum.df.summ, Name) %>%
     arrange(Q.value.Integrated) %>% filter(row_number() <= 1) %>%
     ungroup() %>% arrange(Q.value.Integrated)
write.table(corum.df.summ.unq,'../data/tables/corum.uniqueName.all.tsv', 
            quote = F, sep = '\t', row.names = F)

corumtop30 <- fread('../data/tables/corum.uniqueName.all.tsv')[1:30,] %>%
     mutate(Mean.Integrated.ERC = round(Mean.Integrated.ERC,3),
            Mean.Mammal.ERC = round(Mean.Mammal.ERC,3),
            P.value.Integrated = formatC(P.value.Integrated,3),
            P.value.Mammal = formatC(P.value.Mammal,3),
            Q.value.Integrated = formatC(Q.value.Integrated,3),
            Q.value.Mammal = formatC(Q.value.Mammal,3))

write.table(corumtop30, file = '../data/tables/corum.uniqueName.top30.tsv',
            quote = F, sep = '\t', row.names = F)

length(unique(corum.df.summ$Name))
corum.df.summ.forroc <- select(corum.df.summ, c(4,5))
setnames(corum.df.summ.forroc, c('interc','merc'))
setnames(randcplxdfsumm,c('cplxid','dataset','interc','merc'))
#onembgd <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T) %>% filter(!pair %in% corum.u$df$pair)
head(corum.u$df)

#pdf('../figures/aurocs/corumAllinone.pdf', width = 5, height = 5, family = 'FreeSans')
#dev.off()

#runroc.intvsmamm(corum.u$df, onembgd,plot = T, titl = 'CORUM complex PPIs vs Random PPIs')

pdf('../figures/pvals/corumvsrandomPvalsaurocs.pdf',width = 10,height = 5,family = 'FreeSans')
par(mfrow = c(1,2))
hist(corum.df.summ$P.value.Integrated, breaks = 25, probability = T, xlab = 'p-value', col = rgb(1,0,0,0.3), main = '',
     cex.lab=1.25)
hist(corum.df.summ$P.value.Mammal, breaks = 25, probability = T, col = rgb(0,0,0,0.3), add = T)
legend('right',c('Integrated','Mammal'),title = 'ERC',col = c('red','black'),lty = 1,bty='n',cex=1.25)
title('A',adj = 0, font.main = 1,cex.main = 1.5)
#runroc.intvsmamm(corum.df.summ.forroc, filter(randcplxdfsumm, dataset <= 10),plot = T, titl = 'B')
runpr.intvsmamm(corum.df.summ.forroc, randcplxdfsumm,plot = T, titl = 'B',ypts = c(1,5,10,100,500),ymax=750)
dev.off()

#runroc.intvsmamm(corum.df.summ.forroc, filter(randcplxdfsumm, dataset <= 10000),plot = T, titl = 'CORUM vs Random Complexes')

sum(corum.df.summ$Q.value.Integrated < 0.05, na.rm = T)
sum(corum.df.summ$Q.value.Mammal < 0.05, na.rm = T)


#ignore
nrandoms <- 100
randcplxs <- do.call('rbind',lapply(1:nrandoms, function(y){
     print(y)
     createRandomCplx.df(hsap.genes.useforsim, corumcplx.u$cplx) %>% mutate(dataset = y)
}))
randcplxs$pair <- with(randcplxs,calcpairid(V1,V2,pairmap))
randcplxs <- filter(randcplxs,!pair %in% badblastpairs$pair)
randcplxs$sumlogpfterc <- with(randcplxs,calcpairid(V1,V2,alllogpfterc))
randcplxs$interc <- qnorm(-1*randcplxs$sumlogpfterc, lower.tail = F, log.p = T)
randcplxs$merc <- with(randcplxs,calcpairid(V1,V2,merc.sym))

corumcplx.u$df$lbl <- 1
randcplxs$lbl <- 0

nball <- runNaiveBayes.intvsmamm(corumcplx.u$df, randcplxs,plot=T,titl = 'CORUM Complex PPIs')

corumcplx.u.summ <- group_by(corumcplx.u$df, cplxid) %>%
     summarise(medianscore = mean(interc,na.rm = T))
randcplxs.summ <- group_by(randcplxs, cplxid,dataset) %>%
     summarise(medianscore = mean(interc, na.rm = T))
corumcplx.u.summ$complex <- 'CORUM'
randcplxs.summ$complex <- 'Random'
plotdf <- bind_rows(corumcplx.u.summ, dplyr::select(randcplxs.summ,c(1,3,4)))
plotdf$complex <- factor(plotdf$complex, levels = c('Random','CORUM'))

g <- ggplot(plotdf, aes(x = complex, y=medianscore, fill = complex))+
     geom_violin()+geom_boxplot(width=0.2)+theme_bw()+
     xlab('')+ylab('Mean Integrated Score')+ggtitle('CORUM complex mean score')+
     theme(axis.text=element_text(size = 18, face='bold'),
           axis.title=element_text(size = 20, face='bold'),
           plot.title=element_text(size = 22, face='bold'),
           strip.text = element_text(size = 18, face = 'bold'))
print(g)

corumpvs <- sapply(1:nrow(corumcplx.u.summ), function(x){
     mean(filter(corumcplx.u.summ, cplxid == x)$medianscore<=filter(randcplxs.summ,cplxid == x)$medianscore,
          na.rm = T)
})
quantile(corumpvs, probs = seq(0,1,0.1), na.rm = T)
plot(corumcplx.nunits[corumcplx.nunits>1],-log10(corumpvs))

boxplot(corumcplx.nunits[corumpvs<0.06],corumcplx.nunits[corumpvs>=0.06], ylim = c(0,40))




# nball <- lapply(1:nrandoms, function(y){
#      print(y)
#      runNaiveBayes(corum.erc.df, randcplxs[[y]],plot=F)
# })
# nballfull <- nball
# nball <- nballfull[1:5]
# 
# saveRDS(nballfull, file = '../data/corum.nb.randoms.rds')
# nrandoms <- length(nball)
# nbcorummmdat <- mmdata(scores = join_scores(lapply(nball, '[[', 1), lapply(nball, '[[', 2)), 
#                       labels = join_labels(rep(lapply(nball, '[[', 3),2)), 
#                       modnames = rep(c('Integrated','Mammal only'), each = nrandoms),
#                       dsids = rep(c(1:nrandoms),2))
# mmcurves <- evalmod(nbcorummmdat,cb_alpha = 0.05)
# autoplot(mmcurves, "PRC", show_cb = TRUE )+coord_cartesian(xlim = c(0,0.25))+
#      labs(fill="Method")+
#      theme_bw()+ggtitle('Naive Bayes Prediction: Corum Complexes')+
#      theme(axis.text=element_text(size = 18, face='bold'),
#            axis.title=element_text(size = 20, face='bold'),
#            plot.title=element_text(size = 22, face='bold'),
#            strip.text = element_text(size = 18, face = 'bold'))
