library(RERconverge)

#source R files from RERconverge
source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
source('~/Documents/RERc/RERconverge/R/projection_coevo.R')
source('~/Documents/RERc/RERconverge/R/RcppExports.R')

#source R files for ERC calculations
source('~/Documents/erc/code/funcsERCfromProjections.R')


tres <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
wts <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.weights.rds')

merc <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.erc.rds')
head(sort(merc['NGLY1',], decreasing = T))

nglyv1 <- merc['NGLY1',]
nglyv2 <- merc[,'NGLY1']

nge <- c(nglyv2,nglyv1)
nge <- nge[nge <= 1]
ngeout <- sort(nge, decreasing = T)[3:202]

ngeoutdf <- data.frame(erc = ngeout, gene = names(ngeout), row.names = NULL, stringsAsFactors = F)
write.table(ngeoutdf, file = 'NGLY1.top.erc.gene', row.names = F, col.names = F, quote = F, sep = '\t')

names(head(sort(nge, decreasing = T), 52))

write(head(sort(nge, decreasing = T), 52), file ='ngly1.top.erc')

testset <- c('ZFP90','ZNF189')
testset <- c('SEC61B','GAGE10')
testset <- c('DNAH10','WDR65')
correlateERCGeneListTrees(treesObj = tres, genelist = testset, min.sp = 15, weighted = T, weights = wts)

#calculate relative rates for each gene; row - gene, col - branch
#rers <- getAllResiduals(dtreall,min.sp = 5)
rers <- readRDS('~/Documents/rermethods/data/mamm63nt.trees.rers.sqrt.wt.rds')
plotRers(rers,index = 1) #index can be a number or a gene name

plotRers(rers,index = 'DNAH10') #index can be a number or a gene name
plotRers(rers,index = 'WDR65') #index can be a number or a gene name

#calculate erc of all combinations of pairs of genes; will take a really long time!
ercall <- correlateERCAllTrees(treesObj = dtreall, min.sp = 5)

#erc for a subset pairs of genes 
#define list of genes
dtreall <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')
testset <- names(dtreall$trees)[sample(c(5:10),5)]
erctestset <- correlateERCGeneListTrees(treesObj = tres, genelist = testset, weighted = T, min.sp = 5)
#calculate erc of all combinations of pairs of genes in list
erctestset <- correlateERCGeneListTrees(treesObj = dtreall, genelist = testset, min.sp = 5)
getERCasmat(ercall,testset) #load ERC values for genelist from the pre-computed matrix

par(mfrow = c(1,2))
plot.phylo(tres$trees[['CACNA1B']])
plot.phylo(tres$trees[['CACNA1G']])
