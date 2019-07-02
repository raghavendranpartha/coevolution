nbgd <- 1100000
nbgdf <- 1000000
onembgdpairs <- data.frame(g1 = sample(hsap.genes.useforsim,nbgd,replace = T),
                           g2 = sample(hsap.genes.useforsim,nbgd,replace = T), stringsAsFactors = F) %>%
     mutate(pair = calcpairid(g1,g2,pairmap)) %>%
     filter(pair != 0) %>% 
     mutate(interc = calcpairid(g1,g2,intfterc)) %>%
     mutate(merc = calcpairid(g1,g2,merc.sym)) %>%
     group_by(pair) %>%
     summarise(g1 =g1[1],g2=g2[1],interc=interc[1],merc=merc[1])

sum(is.na(onembgdpairs$interc))
sum(is.na(onembgdpairs$merc))

write.table(onembgdpairs, file = '../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv',
            quote = F, sep = '\t', row.names = F)


fds <- readRDS('../data/fulldatasetintegration_correct.rds')
fds$genepresenceeach[,1:5]

head(10000*fds$genepresenceeach[1,]+
          1000*fds$genepresenceeach[2,]+
          100*fds$genepresenceeach[3,]+
          10*fds$genepresenceeach[4,]+
          fds$genepresenceeach[5,])
genepresenceint <- 10000*fds$genepresenceeach[1,]+1000*fds$genepresenceeach[2,]+100*fds$genepresenceeach[3,]+10*fds$genepresenceeach[4,]+fds$genepresenceeach[5,]

saveRDS(genepresenceint, file = '../data/genepresenceasinteger.rds')

ndatafterc[upper.tri(ndatafterc)] <- NA
diag(ndatafterc) <- NA
ndsmelt <- melt(ndatafterc) %>% na.omit() %>%
     mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))

bgdbig <- rbind(
     sample_n(filter(ndsmelt, value == 1), 5000000),
     sample_n(filter(ndsmelt, value == 2), 5000000),
     sample_n(filter(ndsmelt, value == 3), 5000000)
)
sum(ndsmelt$value == 4)
sum(ndsmelt$value == 5)
bgdbig <- rbind(bgdbig,
                filter(ndsmelt, value == 4),
                filter(ndsmelt, value == 5))
table(bgdbig$value)
rm('ndsmelt')
saveRDS(bgdbig, file = '../data/toperc/bgdbig.rds')

bgdbig <- readRDS('../data/toperc/bgdbig.rds')
pairmap <- readRDS('../data/pairmapmat.rds')
bgdbig$pair <- calcpairid(bgdbig$Var1, bgdbig$Var2, pairmap)
rm(pairmap)
intfterc <- readRDS('../data/allsum.fterc.rds')
bgdbig$interc <- calcpairid(bgdbig$Var1, bgdbig$Var2, intfterc)
rm(intfterc)
merc.sym <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.fterc.symm.rds')
bgdbig$merc <- calcpairid(bgdbig$Var1, bgdbig$Var2, merc.sym)
rm(merc.sym)
     
bgdbig <- filter(bgdbig, !Var1 %in% orgzfg, !Var2 %in% orgzfg)
table(bgdbig$value)
setnames(bgdbig, c('g1','g2','ndataset','pair','interc','merc'))
string.type.use <- fread('../data/public/string.protein.actions.cleaned.tsv', header = T)
bgdbig <- filter(bgdbig, !pair %in% string.type.use$pair)
saveRDS(bgdbig, file = '../data/toperc/bgdbig.rds')
