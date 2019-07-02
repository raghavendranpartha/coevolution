require(dplyr)
require(data.table)
require(ggplot2)

dat <- fread('../data/Marcotte_981complexes.tsv', sep='\t')
#vert.trees <- readRDS('~/Documents/erc/data/saves/vert39.trees.rds')
#write(names(vert.trees$trees), file = '../data/hsap.vert.genes')

cplx <- strsplit(dat$GeneName,';')

allens <- unlist(strsplit(dat$EnsemblID,';'))
datallgenes <- unlist(cplx)

ensidsym <- data.frame(ensid = allens, sym = datallgenes) %>% unique()
length(unique(allens))

hsap.genes <- readLines('../data/hsap.genes')
hsap.vert.genes <- readLines('../data/hsap.vert.genes')
hsap.dmel.genes <- readLines('../data/hsap.dmel.genes')
hsap.cele.genes <- readLines('../data/hsap.cele.genes')
hsap.scer.genes <- readLines('../data/hsap.scer.genes')

cplx.hsap <- sapply(cplx, function(x){
     intersect(x,hsap.genes)
})
cplx.vert <- sapply(cplx, function(x){
     intersect(x,hsap.vert.genes)
})
cplx.dmel <- sapply(cplx, function(x){
     intersect(x,hsap.dmel.genes)
})
cplx.cele <- sapply(cplx, function(x){
     intersect(x,hsap.cele.genes)
})
cplx.scer <- sapply(cplx, function(x){
     intersect(x,hsap.scer.genes)
})

cplx.scer.nunits <- sapply(cplx.scer, length)
cplx.cele.nunits <- sapply(cplx.cele, length)
cplx.dmel.nunits <- sapply(cplx.dmel, length)
cplx.vert.nunits <- sapply(cplx.vert, length)
cplx.hsap.nunits <- sapply(cplx.hsap, length)

table(dat$`Number of Subunits`)
table(cplx.hsap.nunits)

geneomaage <- fread('../data/Marcotte_981complexes.genes.omaage.tsv', header = T) #%>%
     #rowwise() %>%
     #mutate(dist2 = ifelse(distOMA < 0.45,'new','old'))
#dist2 is not useful; checked the presence of unicellular clades in few genes; ProteinAge is the correct prediction

# genesinpaper <- c('TP53','TTN','CSMD3','FAT3','USH2A','MUC16','RYR2','DST','HMCN1','DNAH5','LRP1B','LRP2','APOB','MUC17','ODZ1','RYR1','SYNE1','LRRK2','MACF1','NF1','DNAH3','AHNAK2','PKHD1','SI','DNAH11','GPR98','ANK2','APC','COL22A1','COL6A3','PKHD1L1','PPP1R3A','TG','CRB1','GLI2','RP1L1','FREM2','IGSF10','LAMA2','LAMA3','MAP2','PRKDC','PTPRZ1','STAB2','SYNE2','VPS13B','XIRP2','YSK4','ZFYVE26')

#View(filter(geneomaage, symbol %in% genesinpaper))

marc.cplx.hsap <- unique(unlist(cplx.hsap, use.names = F))     

cplx.age <- sapply(1:length(cplx.hsap), function(x){
     mean(filter(geneomaage, symbol %in% cplx.hsap[[x]])$ProteinAge == 'new')
})
# cplx.age.dist2 <- sapply(1:length(cplx.hsap), function(x){
#      mean(filter(geneomaage, symbol %in% cplx.hsap[[x]])$dist2 == 'new')
# })


cplx.coverage <- data.frame(id = dat$ComplexID,
                            hsap = cplx.hsap.nunits, vert = cplx.vert.nunits, dmel = cplx.dmel.nunits, scer = cplx.scer.nunits, cele = cplx.cele.nunits) %>%
     mutate(vert.pc = vert*100.0/hsap) %>%
     mutate(dmel.pc = dmel*100.0/hsap) %>% mutate(cele.pc = cele*100.0/hsap) %>% 
     mutate(scer.pc = scer*100.0/hsap) %>%
     mutate(avg.pc = (vert.pc+dmel.pc+scer.pc+cele.pc)/4) %>%
     mutate(age = cplx.age) #%>%
     #mutate(age.dist2 = cplx.age.dist2)

#with(cplx.coverage, plot(age, age.dist2))

# with(filter(cplx.coverage, hsap >= 1), plot_colorByDensity(age.dist2, avg.pc))
# abline(lm(avg.pc ~ age, data = cplx.coverage))
with(filter(cplx.coverage, hsap >= 1), plot_colorByDensity(age, avg.pc))
abline(lm(avg.pc ~ age, data = cplx.coverage))

modlm <- lm(avg.pc ~ age, data = cplx.coverage)
summary(modlm)
# modlm2 <- lm(avg.pc ~ age.dist2, data = cplx.coverage)
# summary(modlm2)

g <- ggplot(filter(cplx.coverage, hsap >= 1), aes(age, avg.pc))+
     geom_point()+geom_smooth(method = 'lm')+
     xlab('Complex age (% new proteins)')+ylab('Average coverage\nvert,dmel,scer')+
     theme_bw()+theme(axis.text=element_text(size = 18, face='bold'),
                      axis.title=element_text(size = 20, face='bold'),
                      plot.title=element_text(size = 22, face='bold'),
                      strip.text = element_text(size = 18, face = 'bold'))
print(g)

saveRDS(list(cplx.stats = cplx.coverage,
             cplx=list(hsap = cplx.hsap, vert = cplx.vert, dmel = cplx.dmel, scer = cplx.scer)),
        file = '../data/marcotte981complexes.rds')


plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="",ylab="",main="") {
     
     df <- data.frame(x1,x2)
     x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
     df$dens <- col2rgb(x)[1,] + 1L
     cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
     df$col <- cols[df$dens]
     plot(x2~x1, data=df[order(df$dens),], 
          ylim=ylim,xlim=xlim,pch=20,col=col,
          cex=2,xlab=xlab,ylab=ylab,
          main=main)
}

forbpl <- t(cplx.coverage[,c(6,7,8)])
colnames(forbpl) <- c(1:981)
rownames(forbpl) <- c('vert','dmel','scer')
barplot(forbpl[,c(1:20)], col = c('red','blue','green'), beside = T, xlab = 'Complex ID', ylab = 'Percent genes relative to hsap', legend = rownames(forbpl))
barplot(forbpl[,c(21:40)], col = c('red','blue','green'), beside = T, xlab = 'Complex ID', ylab = 'Percent genes relative to hsap',legend = rownames(forbpl))
barplot(forbpl[,c(41:60)], col = c('red','blue','green'), beside = T, xlab = 'Complex ID', ylab = 'Percent genes relative to hsap',legend = rownames(forbpl))

#cplx.nunits <- sapply(cplx, length)

#plot(mamtres$trees[['RPL17-C18orf32']])
