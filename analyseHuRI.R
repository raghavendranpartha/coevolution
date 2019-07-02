require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)
source('~/Documents/erc/code/funcsERCfromProjections.R')

onembgd <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv')

huri <- fread('../data/public/HuRI.tsv', header = T) %>%
     na.omit()
huri$evi <- rowSums(huri[,c(-1,-2)])

hsap.genes <- readLines('../data/hsap.genes')
ens.sym <- fread('../data/humanSymbolEnsembl.raw.map') %>%
     setnames(c('enstid1','ensgene','enstid2','symbol'))
ens.sym.use <- filter(ens.sym, symbol %in% hsap.genes) %>%
     dplyr::select(c(2,4)) %>% unique()
ens.sym.use.vec <- ens.sym.use$symbol
names(ens.sym.use.vec) <- ens.sym.use$ensgene

huri$genea <- ens.sym.use.vec[huri$Ensembl_gene_id_a]
huri$geneb <- ens.sym.use.vec[huri$Ensembl_gene_id_b]
huri <- dplyr::select(huri, c(1,2,16,17,15,3:14)) %>% na.omit() %>% 
     filter(genea %in% hsap.genes,geneb %in% hsap.genes) %>%
     mutate(pair = calcpairid(genea,geneb,pairmap)) %>%
     filter(pair != 0) %>% group_by(pair) %>%
     filter(row_number()<=1) %>% ungroup()

huri$sumlogpfterc <- calcpairid(huri$genea,huri$geneb,alllogpfterc)
huri$interc <- qnorm(-1*huri$sumlogpfterc,lower.tail = F, log.p = T) 

plotbox2(onembgd$interc, huri$evi, huri$interc, breaksi = c(3,4,6,12), titl = 'HuRI', ylab = 'Integrated score', xlab = 'HuRI evidence score')

huri.erc.df <- spread(do.call('rbind',huri.erc), value = 'value', key = 'dataset')
huri.erc.df$lbl = 1
set.seed(3)
bgdgenelist <- sample(rownames(merc),1000)
bgd.erc.df <- spread(getERCasdf(bgdgenelist), value = 'value', key = 'dataset') %>% mutate(lbl = 0)

library(e1071)
huri.nb.df <- rbind(dplyr::select(huri.erc.df,c(7,2:6)), dplyr::select(bgd.erc.df,c(7,2:6)))
huri.nb.df$lbl <- factor(huri.nb.df$lbl)
Naive_Bayes_Model=naiveBayes(lbl ~., data=huri.nb.df)
Naive_Bayes_Model_mammonly=naiveBayes(lbl ~merc, data=huri.nb.df)
NB_Predictions=predict(Naive_Bayes_Model,huri.nb.df[,c(-1)],type='raw')
NB_Predictions_mammonly=predict(Naive_Bayes_Model_mammonly,huri.nb.df[,c(3)],type='raw')
huri.nb.df$nbp <- NB_Predictions[,2]
huri.nb.df$nbp_mammonly <- NB_Predictions_mammonly[,2]

with(filter(huri.nb.df,lbl == 1), quantile(nbp, probs = seq(0,1,length.out = 10)))
with(filter(huri.nb.df,lbl == 0), quantile(nbp, probs = seq(0,1,length.out = 10)))

require(PRROC)
fg <- filter(huri.nb.df, lbl == 1)$nbp
bg <- filter(huri.nb.df, lbl == 0)$nbp

fg_mammonly <- filter(huri.nb.df, lbl == 1)$nbp_mammonly
bg_mammonly <- filter(huri.nb.df, lbl == 0)$nbp_mammonly


# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr2 <- pr.curve(scores.class0 = fg_mammonly, scores.class1 = bg_mammonly, curve = T)
plot(pr$curve[,1],pr$curve[,2],xlim = c(0,0.1), type = 'l',col='red', lwd = 2, xlab = 'Recall',ylab='Precision',main = 'Naive Bayes Prediction\nhuri PPIs')
lines(pr2$curve[,1],pr2$curve[,2],xlim = c(0,0.1), lwd = 2)
legend('topright', c('Integrated','Mammal only'),lwd = 2,col = c('red','black'),bty='n')

