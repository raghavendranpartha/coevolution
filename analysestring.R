require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)

# string.hg <- fread('../data/public/stringhuman.9606.protein.links.v11.0.txt', header = T) %>%
#      mutate(protein1 = gsub('9606\\.','',protein1)) %>%
#      mutate(protein2 = gsub('9606\\.','',protein2))
string.hg.alias.full <- fread('../data/public/stringhuman9606.protein.aliases.v11.0.txt', header = F)
string.hg.alias <-string.hg.alias.full %>% 
     mutate(V1 = gsub('9606\\.','',V1)) %>%
     filter(grepl('Ensembl_gene',V3)) %>%
     dplyr::select(c(1,2))
# with(string.hg.alias, length(unique(V1)))
# with(string.hg.alias, length(unique(V2)))

hg.ens.sym <- fread('../data/humanSymbolEnsembl.raw.map', header = T) %>%
     dplyr::select(c(2,4)) %>% unique()
string.hg.alias.sym <- left_join(string.hg.alias, hg.ens.sym, by = c('V2'='hg19.ensGene.name2')) %>%
     na.omit() %>% setnames(c('V1','V2','V3')) %>%
     filter(V3 %in% hsap.genes) 

# with(string.hg.alias.sym,length(unique(V1)))
# with(string.hg.alias.sym,length(unique(V3)))
# ensp.sym <- string.hg.alias.sym$V3
# names(ensp.sym) <- string.hg.alias.sym$V1
# head(sort(table(ensp.sym), decreasing = T))
# head(sort(table(string.hg.alias$V3), decreasing = T), 30)
rm('string.hg.alias')
rm('hg.ens.sym')
rm('string.hg.alias.full')
#string.hg.proteins <- with(string.hg,unique(c(protein1,protein2)))
# head(string.hg)
# ens.map <- function(x){
#      ensp.sym[x]
# }
# setDT(string.hg)
# setDT(string.hg.alias.sym)
# string.hg[string.hg.alias.sym,on=.(protein1=V1),g1:=i.V3]
# string.hg[string.hg.alias.sym,on=.(protein2=V1),g2:=i.V3]
# string.hg <- string.hg[!is.na(g1) & !is.na(g2)]
# #string.hg <- string.hg[!is.na(g2)]
# #rm('string.hg')
# #rm('string.hg.use')
# setkey(string.hg,g1)
# string.hg <- string.hg[.(intersect(ensp.sym,hsap.genes))]
# setkey(string.hg,g2)
# string.hg <- string.hg[.(intersect(ensp.sym,hsap.genes))]

# length(intersect(ensp.sym,hsap.genes))
# sum(string.hg$g1 %in% hsap.genes)
# sum(string.hg$g2 %in% hsap.genes)

# rowinds <- match(string.hg$g1, rownames(alllogpfterc))
# colinds <- match(string.hg$g2, colnames(alllogpfterc))
# string.hg$sumlogpfterc <- alllogpfterc[rowinds + nrow(alllogpfterc) * (colinds - 1)]
# 
# write.table(string.hg, file = '../data/public/stringhuman.sumlogpfterc.tsv',
#             quote = F, sep = '\t', row.names = F)

string.type <- fread('../data/public/9606.protein.actions.v11.0.txt', header = T) %>%
     mutate(protein1 = gsub('9606\\.','',item_id_a)) %>%
     mutate(protein2 = gsub('9606\\.','',item_id_b)) %>%
     select(c(8,9,3:7))
setDT(string.type)
setDT(string.hg.alias.sym)
string.type[string.hg.alias.sym,on=.(protein1=V1),g1:=i.V3]
string.type[string.hg.alias.sym,on=.(protein2=V1),g2:=i.V3]
string.type <- string.type[!is.na(g1) & !is.na(g2)]

string.type$interc <- calcpairid(string.type$g1,string.type$g2,intfterc)
string.type$merc <- calcpairid(string.type$g1, string.type$g2, merc.sym)
string.type$pair <- calcpairid(string.type$g1, string.type$g2, pairmap)
string.type$nds <- calcpairid(string.type$g1, string.type$g2, ndatafterc)

string.type.use <- filter(string.type, !g1 %in% orgzfg, !g2 %in% orgzfg,!pair %in% badblastpairs$pair) %>%
     select(c('pair','interc','merc','mode','score','nds')) %>% filter(pair!=0)%>%
     group_by(pair,mode) %>% 
     summarise(interc=interc[1],merc=merc[1],score=max(score),nds=nds[1])

table(string.type$mode)

# toppairs <- fread('../data/toperc/sumlogpftercgt5.bestbitscorepclt10.tsv', header = T)
# toppairs.string <- left_join(toppairs, 
#                              select(string.type, c(3,7,8,9,10,11)),
#                              by = 'pair') %>%
#      filter(!is.na(mode))
# table(toppairs.string$mode)
# length(unique(string.type$pair))
# 
# par(mfrow = c(1,2))
# for(modei in c('activation','binding')){
#      plotbox(filter(toppairs.string, mode == modei)$value, 
#              yy = filter(toppairs.string, mode == modei)$score, 
#              quantiles = F,breaksi = c(seq(5,10,2),25), titl = modei
#      )
# }
# string.type$interclv <- cut(string.type$interc, breaks = c(seq(-6,6,3),25))
# table(filter(string.type, score >= 990)$interclv)
# xx <- table(filter(string.type, score >= 950, mode == 'binding')$interclv)
# xx/sum(xx)
# xx <- table(filter(string.type, score >= 950, mode == 'activation')$interclv)
# xx/sum(xx)
# g <- ggplot(filter(string.type, mode %in% c('activation','binding'),!is.na(interc), score >= 900), aes(x=cut(interc,breaks = c(seq(-6,6,6),25)),y=score))+facet_grid(mode~.)+ geom_violin()
# print(g)

string.type.use.forall <- select(string.type.use, c('pair','interc','score','merc','nds')) %>%
     group_by(pair) %>% summarise(interc=interc[1],merc=merc[1],score=max(score),nds=nds[1])
string.type.use.forall <- string.type.use.forall[!is.na(string.type.use.forall$interc),]
head(sort(table(string.type.use$pair), decreasing = T))
#View(filter(string.type.use.forall, pair %in% c(71513928)))
View(filter(string.type.use, pair %in% c(40845493)))

# onembgdpairs <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T) %>%
#      filter(!pair %in% string.type.use.forall$pair) %>% na.omit()

samplestring.type.use.forall <- group_by(string.type.use.forall, nds) %>%
     filter(row_number() <= 5)
table(string.type.use.forall$nds)

onembgdpairs <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv', header = T) %>% 
     filter(!pair %in% string.type.use.forall$pair)


print(Sys.time())
stringcontrolpairs <- nds.matched.controlpairs(string.type.use.forall$nds, nPerLevelmult = 5)

for(ii in unique(string.type.use$mode)){
     print(ii)
     assign(paste0('stringcontrolpairs.',ii),
            nds.matched.controlpairs(filter(string.type.use, mode == ii)$nds ,nPerLevelmult = 5))
}
print(Sys.time())

plotbox2(stringcontrolpairs,xx = string.type.use.forall$score, 
         yy = string.type.use.forall$interc, nbreaks = 3, titl = paste0('STRING: Any'), 
         xlab = 'Confidence score', ylab = 'Integrated ERC', ylimi=c(-4,11))

pdf('../figures/stringboxplots/any.interc.boxplot.pdf', height = 6, width = 6,family='FreeSans')
plotbox2(onembgdpairs$interc,xx = string.type.use.forall$score, 
         yy = string.type.use.forall$interc, nbreaks = 3, titl = paste0('STRING: Any'), 
         xlab = 'Confidence score', ylab = 'Integrated ERC', ylimi=c(-4,11))
dev.off()
for(i in 1:7){
     ii <- names(table(string.type.use$mode))[i]
     pdf(paste0('../figures/stringboxplots/',ii,'.interc.boxplot.pdf'), 
         height = 6, width = 6,family='FreeSans')
     plotbox2(onembgdpairs$interc,xx = filter(string.type.use, mode==ii)$score, 
              yy = filter(string.type.use, mode==ii)$interc, nbreaks = 3, 
              titl = paste0('STRING: ',ii), xlab = 'Confidence score', ylab = 'Integrated ERC',ylimi=c(-4,11))
     dev.off()
}

pdf('../figures/stringboxplots/any.merc.boxplot.pdf', height = 6, width = 6, family = 'FreeSans')
plotbox2(onembgdpairs$merc,xx = string.type.use.forall$score, 
         yy = string.type.use.forall$merc, nbreaks = 3,ylimi = c(-3,8),
         titl = paste0('STRING: Any'), xlab = 'Confidence score', ylab = 'Mammalian ERC')
dev.off()
for(i in 1:7){
     ii <- names(table(string.type.use$mode))[i]
     pdf(paste0('../figures/stringboxplots/',ii,'.merc.boxplot.pdf'), 
         height = 6, width = 6,family='FreeSans')
     plotbox2(onembgdpairs$merc,xx = filter(string.type.use, mode==ii)$score, 
              yy = filter(string.type.use, mode==ii)$merc, nbreaks = 3, ylimi = c(-3,8),
              titl = paste0('STRING: ',ii), xlab = 'Confidence score', ylab = 'Mammalian ERC')
     dev.off()
}

pdf('../figures/aurocs/stringany.pdf',width = 6,height = 6, family = 'FreeSans')
stringanyrocs <- runroc.string.intvsmamm(true.erc.df = string.type.use.forall, onembgdpairs, levels = c(930,951,960),plot = T ,'STRING: Any',
                                         nPointsforcurve = 1000)
dev.off()
int.levels <- rbind(c(911,940,952),
                    c(943,957,963),
                    c(934,951,958),
                    c(939,952,960))
for(nn in c(1:4)){
     modei = c('activation','binding','catalysis','reaction')[nn]
     pdf(paste0('../figures/aurocs/string',modei,'.pdf'),width = 6,height = 6, family = 'FreeSans')
     assign(paste0(modei,'rocs'),runroc.string.intvsmamm(true.erc.df = as.data.frame(filter(string.type.use, mode == modei)),
                                              onembgdpairs, levels = int.levels[nn,],plot = T ,paste0('STRING: ',modei),nPointsforcurve=1000))
     dev.off()
}
indsformat <-c(3,6,2,5,1,4)
rbind(round(unlist(stringanyrocs, use.names = F)[indsformat],3),
round(unlist(activationrocs, use.names = F)[indsformat],3),
round(unlist(bindingrocs, use.names = F)[indsformat],3),
round(unlist(catalysisrocs, use.names = F)[indsformat],3),
round(unlist(reactionrocs, use.names = F)[indsformat],3))


# runpr.intvsmamm(true.erc.df = filter(string.type.use.forall, score > 960), onembgdpairs, plot = T ,'STRING: Any', xlimi = c(0,1))
# 
# runpr.intvsmamm(df, onembgdpairs, plot = T ,titl = paste0('STRING: ',modei), xlimi = c(0,0.2))
# modei = 'binding'
# df <- as.data.frame(dplyr::filter(string.type.use, mode == modei, score > 963))
# runNaiveBayes.intvsmamm(df, onembgdpairs, plot = T ,paste0('STRING: ',modei))
# runpr.intvsmamm(df, onembgdpairs, plot = T ,titl = paste0('STRING: ',modei), xlimi = c(0,0.2))
