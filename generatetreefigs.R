
mtrep <- mtre
mtrep$masterTree$edge.length <- rep(1,length(mtrep$masterTree$edge.length))
vtrep <- vtre
vtrep$masterTree$edge.length <- rep(1,length(vtrep$masterTree$edge.length))
dtrep <- dtre
dtrep$masterTree$edge.length <- rep(1,length(dtrep$masterTree$edge.length))

ytrep <- ytre
ytrep$masterTree$edge.length <- rep(1,length(ytrep$masterTree$edge.length))

pdf('../figures/trees/mammal.pdf', height = 8,width = 6,family = 'FreeSans')
par(mar = c(0,0,0,0))
plot.phylo(root(mtrep$masterTree,'Platypus'), font = 1,x.lim=c(0,15))
dev.off()

pdf('../figures/trees/vertebrate.pdf', height = 6,width = 4,family = 'FreeSans')
par(mar = c(1,1,1,0))
plot.phylo(root(vtrep$masterTree,node = 64), font = 1)
dev.off()

dmap <- fread('../data/drosophilaspeciesname.map', header = F)
dvecmap <- dmap$V2
names(dvecmap) <- dmap$V1

pdf('../figures/trees/fly.pdf', height = 5,width = 4,family = 'FreeSans')
par(mar = c(1,1,1,0))
dtrep$masterTree$tip.label <- dvecmap[dtrep$masterTree$tip.label]
plot.phylo(dtrep$masterTree, font = 1)
dev.off()

wmap <- fread('../data/wormspeciesname.map', header = F)
wvecmap <- wmap$V2
names(wvecmap) <- wmap$V1

wtrep <- wtre
wtrep$masterTree$edge.length <- rep(1,length(wtrep$masterTree$edge.length))
pdf('../figures/trees/worm.pdf', height = 5,width = 4,family = 'FreeSans')
par(mar = c(1,1,1,0))
wtrep$masterTree$tip.label <- wvecmap[wtrep$masterTree$tip.label]
plot.phylo(wtrep$masterTree, font = 1, x.lim = c(0,14))
dev.off()

sort(ytrep$masterTree$tip.label)

ymap <- fread('../data/yeastspeciesname.map', header = F)
yvecmap <- ymap$V2
names(yvecmap) <- ymap$V1

pdf('../figures/trees/yeast.pdf', height = 5,width = 4,family = 'FreeSans')
par(mar = c(1,1,1,0))
ytrep$masterTree$tip.label <- yvecmap[ytrep$masterTree$tip.label]
plot.phylo(ytrep$masterTree, font = 1)
dev.off()



g1 <- 'TNC'
g2 <- 'DNAH10'

wts <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.weights.rds')
rers <- readRDS('~/Documents/rermethods/data/mamm63nt.trees.rers.sqrt.wt.rds')

source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
source('~/Documents/RERc/RERconverge/R/projection_coevo.R')
source('~/Documents/RERc/RERconverge/R/RcppExports.R')

#source R files for ERC calculations
source('~/Documents/erc/code/funcsERCfromProjections.R')

par(mfrow = c(2,1), mar = c(0,0,0,0))
plot.phylo(root(mtre$trees[[g1]],'Platypus'), x.lim = c(0,0.5), font = 1, cex = 0.8)
plot.phylo(root(mtre$trees[[g2]],'Platypus'), x.lim = c(0,0.5), font = 1, cex = 0.8)

require(gridExtra)
grid.arrange(plotRers(rers,index = g1, xlims = c(-3,2.5)), plotRers(rers,index = g2, xlims = c(-4,3)))

correlateERCGeneListTrees(mtre, c(g1,g2), weights = wts, xlims = c(-2,2))
merc[g1,g2]/sqrt(123-3)
merc[g2,g1]







