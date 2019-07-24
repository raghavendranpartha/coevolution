require(dplyr)
require(data.table)

marc.cplx <- readRDS('../data/public/marcottecplxs.cleaned.rds')
length(marc.cplx$cplx)
min(sapply(marc.cplx$cplx, length))
max(sapply(marc.cplx$cplx, length))

marc.cplxsq <- as.data.frame(t(combn(1:length(marc.cplx$cplx),2)))
marc.cplxsq <- cbind(marc.cplxsq,calc.betweencomplex.interc.numpairs(marc.cplxsq,marc.cplx$cplx))
setnames(marc.cplxsq,c('complexind1','complexind2','meaninterc','numgenepairs'))
marc.cplxsq <- filter(marc.cplxsq, numgenepairs > 2, !is.na(meaninterc))

write.table(marc.cplxsq, file = '../data/public/betweencomplex/marcotte.betweencomplex.tsv', quote = F, sep = '\t', row.names = F)

corum.cplx <- readRDS('../data/public/corumcplxs.cleaned.rds')
corum.cplxsq <- as.data.frame(t(combn(1:length(corum.cplx$cplx),2)))
corum.cplxsq <- cbind(corum.cplxsq,calc.betweencomplex.interc.numpairs(corum.cplxsq,corum.cplx$cplx))
setnames(corum.cplxsq,c('complexind1','complexind2','meaninterc','numgenepairs'))
corum.cplxsq <- filter(corum.cplxsq, numgenepairs > 2, !is.na(meaninterc))
corum.cplxsq$complex1 <- names(corum.cplx$cplx)[corum.cplxsq$complexind1]
corum.cplxsq$complex2 <- names(corum.cplx$cplx)[corum.cplxsq$complexind2]

write.table(corum.cplxsq, file = '../data/public/betweencomplex/corum.betweencomplex.tsv', quote = F, sep = '\t', row.names = F)

head(names(corum.cplx$cplx))

genes <- hsap.genes.useforsim
genepresenceint <- readRDS('../data/genepresenceasinteger.rds')

corum.randcplxs <- generate.randcplxs(corum.cplx$cplx, hsap.genes.useforsim, genepresenceint, nrandoms = 2)
corum.rand.cplxsq <- corum.cplxsq[,c(1:2)]
corum.randcplxs.interc <- do.call('rbind',lapply(1:2, function(y){
     print(y)
     xv <- calc.betweencomplex.interc(corum.rand.cplxsq,corum.randcplxs[[y]])
     cbind(corum.rand.cplxsq,xv) %>% mutate(dataset = y)
}))


# corum <- fread('../data/tables/corum.uniqueName.all.tsv', header = T)
# sum(corum$Q.value.Integrated < 0.05, na.rm = T)
# sum(corum$Q.value.Mammal < 0.05, na.rm = T)
# marcotte <- fread('../data/tables/marcotte.all.tsv')
# sum(marcotte$Q.value.Integrated < 0.05, na.rm = T)
# sum(marcotte$Q.value.Mammal < 0.05, na.rm = T)
# omim <- fread('../data/tables/omim.all.tsv')
# sum(omim$Q.value.Integrated < 0.05, na.rm = T)
# sum(omim$Q.value.Mammal < 0.05, na.rm = T)

corum.betweencomplex <- fread('../data/public/betweencomplex/corum.betweencomplex.withpqvalues.tsv') %>%
     filter(qv < 0.05) 
corum.betweencomplex$complexpair <- with(corum.betweencomplex, 
                                         sapply(1:nrow(corum.betweencomplex), function(x){
                                              paste(sort(c(complex1[x],complex2[x])), collapse = '_')
                                         }))
corum.betweencomplex.use <- corum.betweencomplex %>% group_by(complexpair) %>%
     arrange(desc(numgenepairs)) %>% filter(row_number() <= 1) %>%
     mutate(complexind1 = as.character(complexind1),
            complexind2 = as.character(complexind2))
corumannot <- fread('../data/public/coreComplexes.txt', header = T) %>% 
     filter(Organism == 'Human') %>%
     dplyr::select(c(1,2,11)) 
setnames(corumannot, c('id','complexname','funcatid'))
rawannotmap <- readLines('../data/public/corum.annotation.map')
annotids <- rawannotmap[seq(2,60,by = 3)]
names(annotids) <- rawannotmap[seq(1,60,by = 3)]
corumcomplexannotations <- sapply(1:nrow(corumannot), function(x){
     xids <- strsplit(corumannot$funcatid[x],';')[[1]]
     paste(unique(annotids[sapply(strsplit(xids,'.',fixed = T),'[[',1)]), collapse = ';')
})
names(corumcomplexannotations) <- as.character(corumannot$id)

corum.betweencomplex.use$funcat1 <- corumcomplexannotations[corum.betweencomplex.use$complexind1]
corum.betweencomplex.use$funcat2 <- corumcomplexannotations[corum.betweencomplex.use$complexind2]

corum.betweencomplex.use <- corum.betweencomplex.use %>% na.omit()


