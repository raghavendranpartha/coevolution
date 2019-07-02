annots <- readRDS('../data/annotations.RDS')

names(annots[['mgi']])
head(annots[['human']][['geneset.names']])
head(annots[['mgi']][['geneset.names']])
head(annots[['mgi']][['geneset.descriptions']])

sapply(annots, function(x){
     print(length(x[[2]]))
})
mgi <- annots[['mgi']]
mgi.ngenes <- sapply(mgi[['genesets']], length)
mgi.ngenes.df <- data.frame(ngenes = mgi.ngenes,
                            phen = names(mgi.ngenes), stringsAsFactors = F)
quantile(mgi.ngenes, probs = seq(0,1,0.1))
#sum(mgi.ngenes>2 & mgi.ngenes<150)

mgi.geneset.useraw <- mgi$genesets[mgi.ngenes > 2 & mgi.ngenes <= 150]
# sum(mgi.ngenes > 100)
# sum(mgi.ngenes > 500)

mgi.geneset.use <- sapply(mgi.geneset.useraw, function(x){
     unique(intersect(x, hsap.genes))
})
mgi.geneset.use.n <- sapply(mgi.geneset.use, length)
mgi <- mgi.geneset.use[mgi.geneset.use.n > 1]

mgi.df <- do.call('rbind',lapply(1:length(mgi), function(x){
     as.data.frame(t(combn(mgi[[x]],2))) %>% mutate(cplxid=x)
}))
mgi.df <- mgi.df %>%
     mutate(sumlogpfterc = calcpairid(V1,V2,alllogpfterc)) %>%
     mutate(merc = calcpairid(V1,V2,merc.sym))