require(dplyr)
require(data.table)

irefhg <- fread('../data/public/irefindex/9606.mitab.01-22-2018.txt', header = T)

head(irefhg[,1:5])
head(filter(irefhg,edgetype!='X'),2)
table(irefhg$edgetype)

irefhg$g1 <- gsub('hgnc:','',
     sapply(strsplit(irefhg$aliasA,'|',fixed = T),'[[',1))
irefhg$g2 <- gsub('hgnc:','',
           sapply(strsplit(irefhg$aliasB,'|',fixed = T),'[[',1))

#View(irefhg[1:5,])
irefhg <- irefhg[(irefhg$g1 %in% hsap.genes)&(irefhg$g2 %in% hsap.genes),]
irefhg.use <- select(irefhg,c(12,13,14,15,53,55,56))
rowinds <- match(irefhg.use$g1, rownames(alllogpfterc))
colinds <- match(irefhg.use$g2, colnames(alllogpfterc))
irefhg.use$sumlogpfterc <- alllogpfterc[rowinds + nrow(alllogpfterc) * (colinds - 1)]
irefhg.use$lprconf <- NA
irefhg.use$npconf <- NA
indstodo <- irefhg.use$confidence != '-'
irefhg.use$lprconf[indstodo] <- as.numeric(gsub('lpr:','',
                           sapply(strsplit(irefhg$confidence[indstodo],'|',fixed = T),'[[',2)))
irefhg.use$npconf[indstodo] <- as.numeric(gsub('np:','',
                                      sapply(strsplit(irefhg$confidence[indstodo],'|',fixed = T),'[[',3)))
quantile(filter(irefhg.use,edgetype=='X')$sumlogpfterc, na.rm=T)
irefhg.use$pair <- with(irefhg.use,calcpairid(g1,g2))
write.table(irefhg.use, file = '../data/public/irefindexhuman.sumlogpfterc.tsv',
            quote = F, sep = '\t', row.names = F)

pairmap <- matrix(0,nrow = 19149,ncol=19149,byrow = T)
pairmap[lower.tri(pairmap)] <- 1:choose(19149,2)
pairmap <- pairmap+t(pairmap)
rownames(pairmap) <- rownames(alllogpfterc)
colnames(pairmap) <- rownames(alllogpfterc)
pairmap[1:5,1:5]
saveRDS(pairmap, file = '../data/pairmapmat.rds')