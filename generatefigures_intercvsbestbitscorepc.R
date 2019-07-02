intfterc <- readRDS('../data/allsum.fterc.rds')
intfterc.raw <- readRDS('../data/RAWallsum.fterc.rds')

intfterc.raw[intfterc.raw == 0] <- NA

sum(intfterc[lower.tri(intfterc)] == intfterc.raw[lower.tri(intfterc.raw)], na.rm = T)
sum(lower.tri(intfterc))

diag(intfterc.raw) <- NA
sum(is.na(intfterc.raw))

sum(intfterc[lower.tri(intfterc)] - intfterc.raw[lower.tri(intfterc.raw)], na.rm = T)

onembgd <- fread('../data/toperc/onembgdpairs.bestbitscorepclt0p1.tsv') %>% select(c('g1','g2','interc'))
onembgd$bestbitscorepc <- 0.000001
blastpairs <- fread('../data/allgeneshg19.pairbitscore.tsv') %>%
     mutate(interc = calcpairid(g1,g2,intfterc.raw)) %>%
     select(c('g1','g2','interc','bestbitscorepc'))

plotdf <- rbind(blastpairs, onembgd)
breaks <- c(0,0.0001,0.01,0.05,0.1,0.5,1)

cutres<-cut(plotdf$bestbitscorepc,breaks = breaks,
            right = F)
nbreaks <- length(breaks)
cutres_tt <- table(cutres)



pdf('../figures/intsumftercvsrelativeblastbitscore.pdf', width = 8, height = 4, family = 'FreeSans')
par(mar = c(7,5,1,1))
boxplot((plotdf$interc)~ cutres, xlab = "Relative Blast Bit score", 
        ylab = "Integrated ERC", outline=F,  log="", las=2, srt = 45, xaxt='n', cex.lab = 1.4, cex.axis = 1.25)
abline(h=median(onembgd$interc, na.rm = T), lty = 2, lwd = 2)
axis(1, at=1:(nbreaks), labels = F)
# text(1:(nbreaks-1), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), srt = 45, adj = 1,
#      labels = levels(cutres), xpd = TRUE, cex = 1.5)
text(1:(nbreaks-1), y =-8, srt = 45, adj = 1,
     labels = c('Control',levels(cutres)[-1]), xpd = TRUE, cex = 1.25)
dev.off()
#text(1:(nbreaks), y =1e-07, srt = 30, adj = 0.5,
#     labels = cutres_tt, xpd = TRUE, cex = 1.0)
#title("Integrated ERC vs Relative Blast bit score")