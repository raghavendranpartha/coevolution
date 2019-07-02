library(RERconverge)
require(dplyr)
require(data.table)
hsap.genes <- readLines('../data/hsap.genes')
hsap.vert.genes <- readLines('../data/hsap.vert.genes')
hsap.dmel.genes <- readLines('../data/hsap.dmel.genes')
hsap.cele.genes <- readLines('../data/hsap.cele.genes')
hsap.scer.genes <- readLines('../data/hsap.scer.genes')

mbeta <- readRDS('../data/ercbetafit/mammalercbetafit.rds')$est
vbeta <- readRDS('../data/ercbetafit/vertercbetafit.rds')$est
wbeta <- readRDS('../data/ercbetafit/wormercbetafit.rds')$est
dbeta <- readRDS('../data/ercbetafit/dmelercbetafit.rds')$est
ybeta <- readRDS('../data/ercbetafit/scerercbetafit.rds')$est


genepresence <- (hsap.genes %in% hsap.vert.genes)+(hsap.genes %in% hsap.dmel.genes)+(hsap.genes %in% hsap.scer.genes)+(hsap.genes %in% hsap.cele.genes)
names(genepresence) = hsap.genes
table(genepresence)
genesinmorethanonedataset <- hsap.genes[which(genepresence == 4)]
genesin3 <- hsap.genes[which(genepresence == 3)]
write(genesinmorethanonedataset, file = '../data/common.all5lineages.gene')
write(genesin3, file = '../data/common.4lineages.gene')
#commonpairs <- as.data.frame(t(combn(genesinmorethanonedataset,2)), stringsAsFactors = F)
commonpairs <- as.data.frame(t(combn(genesin3,2)), stringsAsFactors = F)
commonpairs$gene1c = genepresence[commonpairs$V1]
commonpairs$gene2c = genepresence[commonpairs$V2]
commonpairs$pairc = with(commonpairs, gene1c+gene2c)
commonpairs <- arrange(commonpairs, desc(pairc))

head(commonpairs, 20)

npairs <- 1000
testpairs <- commonpairs[1:npairs,]



#commongenes <- intersect(intersect(hsap.vert.genes,hsap.dmel.genes),hsap.scer.genes)

flymap <- fread('../data/dmel.hsap.geneIDsUniProtIDs.map')
yeastmap <- fread('../data/hsap.scer.geneIDsUniProtIDs.map')
wormmap <- fread('../data/cele.OgwbgnUniprotSymbol.map')
wormmap.vec <- wormmap$og
names(wormmap.vec) <- wormmap$symbol
flymap.vec <- flymap$og
names(flymap.vec) <- flymap$hsapsymbol
yeastmap.vec <- yeastmap$scersymbol
names(yeastmap.vec) <- yeastmap$hsapsymbol

saveRDS(list(flymap = flymap.vec,
             wormmap = wormmap.vec,
             yeastmap = yeastmap.vec), file = '../data/flywormyeastmap.vec.rds')

getercasvec <- function(mat, gen){
     #mat <- merc
     genv <- match(gen,rownames(mat))
     genord <- gen[order(genv)]
     out <- mat[genord,genord][lower.tri(mat[genord,genord], diag = F)]
     out[is.infinite(out)] <- NA
     out
}

commonpairserc <- sapply(c(1:nrow(commonpairs)), function(x){
     #x <- 1
     #gen1 <- commonpairs[x,1]
     #gen2 <- commonpairs[x,2]
     #x <- 1
     c(getercasvec(merc, c(commonpairs$V1[x],commonpairs$V2[x])),
       getercasvec(verc, c(commonpairs$V1[x],commonpairs$V2[x])),
       getercasvec(werc, wormmap.vec[c(commonpairs$V1[x],commonpairs$V2[x])]),
       getercasvec(derc, flymap.vec[c(commonpairs$V1[x],commonpairs$V2[x])]),
       getercasvec(yerc, yeastmap.vec[c(commonpairs$V1[x],commonpairs$V2[x])]))
})

commonpairs.score <- cbind(commonpairs, t(commonpairserc)) %>%
     setnames(c('gene1','gene2','gene1c','gene2c','pairc','merc','verc','werc','derc','yerc')) %>%
     mutate(msc = -log(pbeta(merc,mbeta[1],mbeta[2],lower.tail = F),10)) %>%
     mutate(vsc = -log(pbeta(verc,vbeta[1],vbeta[2],lower.tail = F),10)) %>%
     mutate(wsc = -log(pbeta(werc,wbeta[1],wbeta[2],lower.tail = F),10)) %>%
     mutate(dsc = -log(pbeta(derc,dbeta[1],dbeta[2],lower.tail = F),10)) %>%
     mutate(ysc = -log(pbeta(yerc,ybeta[1],ybeta[2],lower.tail = F),10)) %>%
     rowwise() %>% mutate(sc = sum(c(msc,vsc,wsc,dsc,ysc), na.rm = T)) %>% arrange(desc(sc)) %>%
     select(c(1:5,16,6:15))

write.table(commonpairs.score, file = '../data/common.all5lineages.pairs.erc.tsv',
            quote = F, row.names = F, sep = '\t')

sumerc <- colSums(commonpairserc, na.rm = T)
quantile(sumerc)
sum(sumerc > 1)

saveRDS(list(genes = commongenes, pairs = commonpairs, erc = commonpairserc), file = 'CommonGenePairsErc.rds')