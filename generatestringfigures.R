string.type.use <- fread('../data/public/string.protein.actions.cleaned.tsv', header = T)
string.type.use.forall <- select(string.type.use, c('pair','interc','score','merc','ndataset')) %>%
     group_by(pair) %>% summarise(interc=interc[1],merc=merc[1],score=max(score),ndataset=ndataset[1])

bgdbig <- readRDS('../data/toperc/bgdbig.rds')
set.seed(2)
#stringany.controlpairs <- ndataset.matched.controlpairs(filter(string.type.use.forall,score > 928)$ndataset, nPerLevelmult = 10)
qs.any <- quantile(string.type.use.forall$score, c(0.9,0.95,0.98,0.99,1))
stringany.controlpairs.with10pc <- ndataset.matched.controlpairs(filter(string.type.use.forall,score > 913)$ndataset, nPerLevelmult = 10)

for(ii in 1:(length(qs.any)-1)){
     print(wilcox.test(filter(string.type.use.forall, score > qs.any[ii], score <= qs.any[ii+1])$interc,
                 stringany.controlpairs.with10pc$interc,
                 alternative = 'g')$p.v)
}
for(ii in 1:(length(qs.any)-1)){
     print(wilcox.test(filter(string.type.use.forall, score > qs.any[ii])$interc,
                       stringany.controlpairs.with10pc$interc,
                       alternative = 'g')$p.v)
}
for(ii in 1:(length(qs.any)-2)){
     print(paste0(qs.any[ii+1], qs.any[ii+2]))
     print(paste0(qs.any[ii], qs.any[ii+1]))
     print(wilcox.test(filter(string.type.use.forall, score > qs.any[ii+1], score <= qs.any[ii+2])$interc,
                       filter(string.type.use.forall, score > qs.any[ii], score <= qs.any[ii+1])$interc,
                       alternative = 'g')$p.value)
}

top10521levels <- sapply(unique(string.type.use$mode), function(xx){
     scs <- filter(string.type.use, mode == xx)$score
     quantile(scs, probs = c(0.9,0.95,0.98,0.99), na.rm = T)
})
rownames(top10521levels) <- c('top10','top5','top2','top1')   

for(ii in unique(string.type.use$mode)){
     print(ii)
     assign(paste0('stringcontrolpairs.',ii),
            ndataset.matched.controlpairs(filter(string.type.use, mode == ii, score > top10521levels['top10',ii])$ndataset ,nPerLevelmult = 10))
}

pdf('../figures/stringboxplots/any.interc.boxplot.pdf', height = 6, width = 8,family='FreeSans')
# plotbox2(stringany.controlpairs$interc,xx = string.type.use.forall$score, 
#          yy = string.type.use.forall$interc, nbreaks = 3, titl = paste0('STRING PPIs'), 
#          xlab = 'Confidence score', ylab = 'Integrated ERC', ylimi=c(-4,11))
plotbox2(stringany.controlpairs.with10pc$interc,xx = string.type.use.forall$score, 
         yy = string.type.use.forall$interc, breaksi = qs.any, leglabsi = c('Control','Top 10-5%','Top 5-2%','Top 2-1%','Top 1%'), titl = paste0('STRING PPIs'), xlab = 'Confidence score', ylab = 'Integrated ERC', ylimi=c(-4,11))
dev.off()
for(i in 1:7){
     ii <- names(table(string.type.use$mode))[i]
     pdf(paste0('../figures/stringboxplots/',ii,'.interc.boxplot.pdf'), 
         height = 6, width = 8,family='FreeSans')
     plotbox2(get(paste0('stringcontrolpairs.',ii))$interc,xx = filter(string.type.use, mode==ii)$score, 
              yy = filter(string.type.use, mode==ii)$interc, breaksi = c(top10521levels[,ii],996),
              leglabsi = c('Control','Top 10-5%','Top 5-2%','Top 2-1%','Top 1%'),
              titl = paste0('Interaction mode: ',ii), xlab = 'Confidence score', ylab = 'Integrated ERC',ylimi=c(-4,11))
     dev.off()
}

pdf('../figures/aurocs/stringany.pdf',width = 6,height = 6, family = 'FreeSans')
# stringanyrocs <- runroc.string.intvsmamm(true.erc.df = string.type.use.forall, stringany.controlpairs, 
#                                          levels = c(928,951,960),plot = T ,'A. Interaction mode: Any',
#                                          nPointsforcurve = 1000)
stringanyrocs <- runroc.string.intvsmamm(true.erc.df = string.type.use.forall, stringany.controlpairs.with10pc, 
                                         levels = c(913,928,951,960),plot = T ,'A. Interaction mode: Any',
                                         nPointsforcurve = 1000)
dev.off()

aucs.intandmamm <- sapply(unique(string.type.use$mode), function(modei){
     runroc.string.intvsmamm(true.erc.df = as.data.frame(filter(string.type.use, mode == modei)),
                             get(paste0('stringcontrolpairs.',modei)), 
                             levels = top10521levels[,modei],plot = F)
})


aucs.intandmamm.forany <- runroc.string.intvsmamm(true.erc.df = as.data.frame(filter(string.type.use.forall, score > 913)),
                             stringany.controlpairs.with10pc, levels = qs.any[-length(qs.any)],plot = F)

dfaucs.intandmamm.forany <- data.frame(auc = unlist(aucs.intandmamm.forany, use.names = F),
                                       Mode = 'Any',
                                       Method = rep(c('A. Integrated ERC','B. Mammal ERC'), each = 4),
                                       Dataset = factor(rep(c('Top10%','Top 5%','Top 2%','Top 1%'), 2), 
                                                        levels = c('Top10%','Top 5%','Top 2%','Top 1%')))

aucdf <- data.frame(auc = unlist(aucs.intandmamm, use.names = F), 
                    Mode = factor(rep(unique(string.type.use$mode), each = 8)),
                    Method = rep(rep(c('A. Integrated ERC','B. Mammal ERC'), each = 4),7),
                    Dataset = factor(rep(c('Top10%','Top 5%','Top 2%','Top 1%'), 14), 
                                     levels = c('Top10%','Top 5%','Top 2%','Top 1%')))


plotdf <- rbind(aucdf, dfaucs.intandmamm.forany)
plotdf$isany <- 1*(plotdf$Mode == 'Any')

require(ggplot2)
pdf('../figures/aurocs/string.modes.auc.top521.intvsmammerc.pdf', width = 10,height = 5,family = 'FreeSans')
g <- ggplot(filter(plotdf,Mode!='catalysis'), 
            aes(x=Dataset,y=auc, shape=Mode, 
                group = Mode, color=factor(isany)))+
     facet_grid(.~Method)+
     scale_color_manual(values=c("black","red"), guide = F)+
     scale_shape_manual(values=c(0:6,9))+
     coord_cartesian(ylim=c(0.45,0.75))+
     ylab('AUC')+geom_point(size=3,stroke = 1.5)+geom_line(linetype='solid')+theme_bw()+
     theme(axis.text=element_text(size = 16,colour = 'black'),
           axis.title=element_text(size = 18),
           plot.title=element_text(size = 22),
           strip.text = element_text(size = 18),
           legend.title=element_text(size = 18),
           legend.text=element_text(size = 16))
print(g)
dev.off()


# int.levels <- rbind(c(911,940,952),
#                     c(942,957,963),
#                     c(932,950,958),
#                     c(936,952,960))
# for(nn in c(1:4)){
#      modei = c('activation','binding','catalysis','reaction')[nn]
#      pdf(paste0('../figures/aurocs/string',modei,'.pdf'),width = 6,height = 6, family = 'FreeSans')
#      assign(paste0(modei,'rocs'),runroc.string.intvsmamm(true.erc.df = as.data.frame(filter(string.type.use, mode == modei)),
#                                                          get(paste0('stringcontrolpairs.',ii)), levels = int.levels[nn,],plot = T ,
#                                                          paste0('B. Interaction mode: ',modei),nPointsforcurve=1000))
#      dev.off()
# }
# indatasetformat <-c(3,6,2,5,1,4)
# rbind(round(unlist(stringanyrocs, use.names = F)[indatasetformat],3),
#       round(unlist(activationrocs, use.names = F)[indatasetformat],3),
#       round(unlist(bindingrocs, use.names = F)[indatasetformat],3),
#       round(unlist(catalysisrocs, use.names = F)[indatasetformat],3),
#       round(unlist(reactionrocs, use.names = F)[indatasetformat],3))
