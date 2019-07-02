require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)
library(RERconverge)

library(fitdistrplus)

marc.cplx <- readRDS('../data/marcotte981complexes.rds')
View(marc.cplx$cplx.stats)
oldcplxinds <- which(marc.cplx$cplx.stats$age == 0)
newcplxinds <- which(marc.cplx$cplx.stats$age > 0)
geneomaage <- fread('../data/Marcotte_981complexes.genes.omaage.tsv', header = T)
#sum(is.na(geneomaage$symbol))

merc <- readRDS('~/Documents/erc/data/saves/mamm63.erc.rds')
#sum(is.infinite(merc))
#sum(is.infinite(verc))
#sum(is.infinite(derc))
#sum(is.infinite(yerc))
verc <- readRDS('~/Documents/erc/data/saves/vert39.erc.rds')
derc <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.erc.rds')
yerc <- readRDS('~/Documents/erc/data/saves/all_Scer.erc.rds')

plotbgderc <- function(mat,npairs = 100000,titl = ''){
     dat <- sample(mat[lower.tri(mat)], npairs)
     dat[is.infinite(dat)] <- NA
     denserc <- density(dat, na.rm = T)
     vs <- quantile(dat, c(0.05,0.95),na.rm = T)
     plot(denserc,main = titl)
     abline(v=vs[1],lty = 2)
     abline(v=vs[2],lty = 2)
}
par(mfrow = c(2,2))
plotbgderc(merc,titl = 'Mammal')
plotbgderc(verc,titl = 'Vert')
plotbgderc(derc,titl = 'Dmel')
plotbgderc(yerc,titl = 'Scer')



set.seed(2)
numpairs <- 500000
dat <- sample(merc[lower.tri(merc)], numpairs)
dat[is.infinite(dat)] <- NA
denserc <- density(dat, na.rm = T)
datscaled <- (dat+1)/2
betafit <- fitdist(datscaled[!is.na(datscaled)], distr = 'beta',method = 'mle',lower = c(0, 0))
predict(betafit,c(0,0.5,0.75))
hist(sample(merc[lower.tri(merc)], numpairs), breaks = 25, probability = T, ylim = c(0,max(denserc[['y']])))
lines(denserc)

flymap <- fread('../data/dmel.hsap.geneIDsUniProtIDs.map')
yeastmap <- fread('../data/hsap.scer.geneIDsUniProtIDs.map')
flymap.vec <- flymap$og
names(flymap.vec) <- flymap$hsapsymbol
yeastmap.vec <- yeastmap$scersymbol
names(yeastmap.vec) <- yeastmap$hsapsymbol

head(rownames(merc))
head(rownames(verc))
head(rownames(derc))
head(rownames(yerc))

getercasvec <- function(mat, gen){
     #mat <- merc
     genv <- match(gen,rownames(mat))
     genord <- gen[order(genv)]
     out <- mat[genord,genord][lower.tri(mat[genord,genord], diag = F)]
     out[is.infinite(out)] <- NA
     out
}

cplx.erc <- lapply(1:nrow(marc.cplx$cplx.stats), function(x){
     #x <- 150
     #gen <- marc.cplx$cplx$hsap[[x]]
     list(merc = getercasvec(merc, marc.cplx$cplx$hsap[[x]]),
          verc = getercasvec(verc, marc.cplx$cplx$vert[[x]]),
          derc = getercasvec(derc, flymap.vec[marc.cplx$cplx$dmel[[x]]]),
          yerc = getercasvec(yerc, yeastmap.vec[marc.cplx$cplx$scer[[x]]]))
     
})

cplx.erc.asvec <- unlist(cplx.erc, use.names = F)
cplx.erc.asvec.hsap <- unlist(sapply(1:length(cplx.erc), function(x){cplx.erc[[x]][[1]]}), use.names = F)
cplx.erc.asvec.vert <- unlist(sapply(1:length(cplx.erc), function(x){cplx.erc[[x]][[2]]}), use.names = F)
cplx.erc.asvec.dmel <- unlist(sapply(1:length(cplx.erc), function(x){cplx.erc[[x]][[3]]}), use.names = F)
cplx.erc.asvec.scer <- unlist(sapply(1:length(cplx.erc), function(x){cplx.erc[[x]][[4]]}), use.names = F)
hist(cplx.erc.asvec, breaks = 25, probability = T)
plot(denserc, lty = 2, col = 'black', 'ERC distribution Marcotte Complex PPIs', xlab = 'ERC')
abline(v=0, lty = 1)
#lines(density(cplx.erc.asvec, na.rm = T))
lines(density(cplx.erc.asvec.hsap, na.rm = T), col = 'red')
lines(density(cplx.erc.asvec.vert, na.rm = T), col = 'blue')
lines(density(cplx.erc.asvec.dmel, na.rm = T), col = 'green')
lines(density(cplx.erc.asvec.scer, na.rm = T), col = 'brown')
legend('topright',c('Mammal background','Mammal','Vert','Dmel','Scer'),col=c('black','red','blue','green','brown'), lty = c(2,1,1,1,1))

oldcplx.erc.asvec.hsap <- unlist(sapply(oldcplxinds, function(x){cplx.erc[[x]][[1]]}), use.names = F)
oldcplx.erc.asvec.vert <- unlist(sapply(oldcplxinds, function(x){cplx.erc[[x]][[2]]}), use.names = F)
oldcplx.erc.asvec.dmel <- unlist(sapply(oldcplxinds, function(x){cplx.erc[[x]][[3]]}), use.names = F)
oldcplx.erc.asvec.scer <- unlist(sapply(oldcplxinds, function(x){cplx.erc[[x]][[4]]}), use.names = F)
newcplx.erc.asvec.hsap <- unlist(sapply(newcplxinds, function(x){cplx.erc[[x]][[1]]}), use.names = F)
newcplx.erc.asvec.vert <- unlist(sapply(newcplxinds, function(x){cplx.erc[[x]][[2]]}), use.names = F)
newcplx.erc.asvec.dmel <- unlist(sapply(newcplxinds, function(x){cplx.erc[[x]][[3]]}), use.names = F)
newcplx.erc.asvec.scer <- unlist(sapply(newcplxinds, function(x){cplx.erc[[x]][[4]]}), use.names = F)

par(mfrow = c(1,2))
plot(denserc, lty = 2, col = 'black', 'ERC distribution Marcotte Old Complexes PPIs', xlab = 'ERC')
abline(v=0, lty = 2)
abline(v=0.25, lty = 2)
abline(v=-0.25, lty = 2)
#lines(density(cplx.erc.asvec, na.rm = T))
lines(density(oldcplx.erc.asvec.hsap, na.rm = T), col = 'red')
lines(density(oldcplx.erc.asvec.vert, na.rm = T), col = 'blue')
lines(density(oldcplx.erc.asvec.dmel, na.rm = T), col = 'green')
lines(density(oldcplx.erc.asvec.scer, na.rm = T), col = 'brown')
# legend('topright',c('Mammal background','Mammal','Vert','Dmel','Scer'),col=c('black','red','blue','green','brown'), lty = c(2,1,1,1,1))
plot(denserc, lty = 2, col = 'black', 'ERC distribution Marcotte New Complexes PPIs', xlab = 'ERC')
abline(v=0, lty = 2)
abline(v=0.25, lty = 2)
abline(v=-0.25, lty = 2)
#lines(density(cplx.erc.asvec, na.rm = T))
lines(density(newcplx.erc.asvec.hsap, na.rm = T), col = 'red')
lines(density(newcplx.erc.asvec.vert, na.rm = T), col = 'blue')
lines(density(newcplx.erc.asvec.dmel, na.rm = T), col = 'green')
lines(density(newcplx.erc.asvec.scer, na.rm = T), col = 'brown')
# legend('topright',c('Mammal background','Mammal','Vert','Dmel','Scer'),col=c('black','red','blue','green','brown'), lty = c(2,1,1,1,1))

cplx.erc.mean <- data.frame (id = c(1:981), stringsAsFactors = F) %>%
     mutate(mamm = sapply(1:length(cplx.erc), function(x) {mean(cplx.erc[[x]][[1]], na.rm = T)})) %>%
     mutate(vert = sapply(1:length(cplx.erc), function(x) {mean(cplx.erc[[x]][[2]], na.rm = T)})) %>%
     mutate(dmel = sapply(1:length(cplx.erc), function(x) {mean(cplx.erc[[x]][[3]], na.rm = T)})) %>%
     mutate(scer = sapply(1:length(cplx.erc), function(x) {mean(cplx.erc[[x]][[4]], na.rm = T)}))  %>%
     rowwise() %>%
     mutate(meanall = mean(c(mamm,vert,dmel,scer), na.rm = T))

hist(cplx.erc)

hist(merc[lower.tri(merc)], breaks = 25, probability = T)
plot(density(merc[lower.tri(merc)], breaks = 25))
     
hist(cplx.erc.mean$mamm, breaks = 25, col = rgb(1,0,0,0.2), probability = T)
hist(cplx.erc.mean$vert, breaks = 25, col = rgb(0,1,0,0.2), add = T, probability = T)
hist(cplx.erc.mean$dmel, breaks = 25, col = rgb(0,0,1,0.2), add = T, probability = T)
hist(cplx.erc.mean$scer, breaks = 25, col = rgb(1,0,1,0.2), add = T, probability = T)


