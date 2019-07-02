require(ggplot2)
string.hg <- fread('../data/public/stringhuman.sumlogpfterc.tsv', sep = '\t', header = T)
# string.hg$pair <- with(string.hg, calcpairid(g1,g2))
# write.table(string.hg, file = '../data/public/stringhuman.sumlogpfterc.tsv',
#              quote = F, sep = '\t', row.names = F)
toppairs <- fread('../data/toperc/sumlogpftercgt5.bestbitscorepclt10.tsv', header = T)
toppairs.string <- inner_join(toppairs, filter(string.hg, pair %in% toppairs$pair), by = 'pair')
#toppairs.string$combined_score[is.na(toppairs.string$combined_score)] <- 0

with(toppairs.string, plotbox(sumlogpfterc,combined_score,quantiles = F,nbreaks = 5,
                              breaksi = c(5,7.5,10,20,75),
                              ylab = 'STRING score',xlab = 'Integrated ERC score'))

toppairs.string$valuelv <- cut(toppairs.string$value, breaks = c(0,15,16,17,19,27,75,100))
g <- ggplot(toppairs.string,aes(x=valuelv,y=combined_score))+
     coord_cartesian(ylim = c(150,1000))+
     geom_violin(adjust=0.1,yScale='log10')+geom_boxplot(width = 0.1)
print(g)


boxplot(filter(string.hg, combined_score >= 980)$sumlogpfterc, prob = T, na.rm = T,log='y')
boxplot(log(filter(string.hg, combined_score >= 998)$sumlogpfterc),
        log(filter(string.hg, combined_score <= 200)$sumlogpfterc),
        prob = T, na.rm = T,log='')
quantile(filter(string.hg, combined_score == 150)$sumlogpfterc, na.rm = T)
quantile(filter(string.hg, combined_score >= 998)$sumlogpfterc, na.rm = T)
#quantile(filter(string.hg, combined_score >= 990)$sumlogpfterc, na.rm = T)

plotbox

stringpairs1 <- with(string.hg,paste(g1,g2,sep=','))
#stringpairs2 <- with(string.hg,paste(g2,g1,sep=','))
head(toppairs$pair %in% stringpairs1)
#head(toppairs$pair %in% stringpairs2)
string.toppairs <- filter(string.h)

string.top <- filter(string.hg, combined_score>960)
rowinds <- match(string.top$g1, rownames(alllogpfterc))
colinds <- match(string.top$g2, colnames(alllogpfterc))
string.top$sumlogpfterc <- alllogpfterc[rowinds + nrow(alllogpfterc) * (colinds - 1)]
string.toppairs1 <- with(string.top,paste(g1,g2,sep=','))
sum(toppairs$pair %in% string.toppairs1)

inds <- commoninds(string.top[,c(4,5)],toppairs[,c(2,3)])
View(string.top[inds,])


toppairs <- fread('../data/toperc/sumlogpftercgt15.bestbitscorepclt10.tsv',header = T)

