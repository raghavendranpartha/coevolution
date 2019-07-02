library(RERconverge)
require(dplyr)
require(data.table)

mamtrees <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
write(names(mamtrees$trees), file = '../data/hsap.genes')
head(names(mamtrees$trees))

hsap.map <- fread('../data/genesymbolUCSCnameUniProt.raw.map', header = T) %>%
     setnames(c('ucscid','protid','kgprotid','dispid','symbol'))


hg19.uniprot.wdmelortholog <- readLines('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens.hsap.uniprot.id')
hg19.uniprot.wscerortholog <- readLines('../data/inparanoid.ortholog.maps/hsap_scer/cleanedparsedOutput.H.sapiens-S.cerevisiae.hsap.uniprot.id')
hg19.uniprot.wceleortholog <- readLines('../data/inparanoid.ortholog.maps/cele_hsap/cleanedparsedOutput.C.elegans-H.sapiens.hsap.uniprot.id')
allhg19uniprot2 <- unique(c(hg19.uniprot.wdmelortholog, hg19.uniprot.wscerortholog))
allhg19uniprot <- unique(c(allhg19uniprot,hg19.uniprot.wceleortholog))

sum(allhg19uniprot2 %in% hsap.map$protid)

with(hsap.map, identical(protid, kgprotid))
noma <- with(hsap.map, which(protid != kgprotid))
View(hsap.map[noma,])
hsap.map.use <- hsap.map %>% dplyr::select(c(3,5)) %>% unique() %>%
     filter(kgprotid %in% allhg19uniprot) %>%
     filter(symbol %in% names(mamtrees$trees)) %>%
     setnames(c('uniprot','symbol'))

write.table(hsap.map.use, file = '../data/hsap.uniprotsymbol.map',
            quote = F, row.names = F, sep = '\t')

sum(hsap.map.use$symbol %in% names(mamtrees$trees))


#map scer
hsap.map.use <- fread('../data/hsap.uniprotsymbol.map', header = T)
mamtrees <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
ytres <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.rds')

hsap.scer.ortho <- fread('../data/inparanoid.ortholog.maps/hsap_scer/cleanedparsedOutput.H.sapiens-S.cerevisiae') %>%
     select(c(2,3)) %>% setnames(c('hsapuniprot','sceruniprot'))

scer.uniprot <- read.fwf('../data/yeast_uniprot.map.txt', widths = c(75,9,11,6,40), header = F, stringsAsFactors = F) %>%
     dplyr::select(c(2,4)) %>%
     mutate(V2 = trim(V2)) %>%
     filter(V2 %in% names(ytres$trees))

#sum(hsap.scer.ortho$sceruniprot %in% scer.fbg.uniprot$V4)
#sum( %in% hsap.scer.ortho$sceruniprot)

hsap.scer.full <- left_join(hsap.scer.ortho, scer.uniprot, by = c('sceruniprot'='V4')) %>%
     left_join(hsap.map.use, by = c('hsapuniprot'='uniprot')) %>% na.omit() %>%
     setnames(c('hsapuniprot','sceruniprot','scersymbol','hsapsymbol'))
sum(hsap.scer.full$scersymbol %in% names(ytres$trees))
sum(hsap.scer.full$hsapsymbol %in% names(mamtrees$trees))

write.table(hsap.scer.full, file = '../data/hsap.scer.geneIDsUniProtIDs.map', quote = F, 
            row.names = F, sep = '\t')

hsap.scer.full <- fread('../data/hsap.scer.geneIDsUniProtIDs.map') 
hsap.scer.full.unqhsapsymbol <- hsap.scer.full %>%
     group_by(hsapsymbol) %>%
     filter(row_number() <= 1) %>%
     ungroup()

write.table(hsap.scer.full.unqhsapsymbol, file = '../data/hsap.scer.geneIDsUniProtIDs.map', quote = F, 
            row.names = F, sep = '\t')

write(hsap.scer.full.unqhsapsymbol$hsapsymbol, file = '../data/hsap.scer.genes')

length(hsap.scer.full$hsapsymbol)
length(unique(hsap.scer.full$hsapsymbol))

head(sort(table(hsap.scer.full$hsapsymbol), decreasing = T))

#map dmel
hsap.map.use <- fread('../data/hsap.uniprotsymbol.map', header = T)
mamtrees <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
dtres <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')

dmel.hsap.ortho <- fread('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens') %>% select(c(2,3)) %>% setnames(c('dmeluniprot','hsapuniprot'))

dmel.fbgsymunip <- fread('../data/dmel.OgFbgnUniprotSymbol.map', header = T)

dmel.hsap.full <- left_join(dmel.hsap.ortho, hsap.map.use, by = c('hsapuniprot'='uniprot')) %>%
     left_join(dmel.fbgsymunip, by = c('dmeluniprot'='uniprot')) %>% na.omit() %>%
     setnames(c('dmeluniprot','hsapuniprot','hsapsymbol','fbgn','og','dmelsymbol'))

write.table(dmel.hsap.full, file = '../data/dmel.hsap.geneIDsUniProtIDs.map', quote = F,
            row.names = F, sep = '\t')
dmel.hsap.full <- fread('../data/dmel.hsap.geneIDsUniProtIDs.map', header = T)
dmel.hsap.full[which(!dmel.hsap.full$og %in% rownames(erc)),] # erc is the droso erc

write.table(dmel.hsap.full[which(!dmel.hsap.full$og %in% rownames(erc)),], file = '../data/dmel.hsap.geneIDsUniProtIDs.map', quote = F,
            row.names = F, sep = '\t')

write(dmel.hsap.full[which(!dmel.hsap.full$og %in% rownames(erc)),]$hsapsymbol, file = '../data/hsap.dmel.genes')

sum(dmel.hsap.full$og %in% names(dtres$trees))
length(dmel.hsap.full$hsapsymbol)
length(unique(dmel.hsap.full$hsapsymbol))

#map cele
hsap.map.use <- fread('../data/hsap.uniprotsymbol.map', header = T)
length(hsap.map.use$symbol)
length(unique(hsap.map.use$symbol))
length(unique(hsap.map.use$uniprot))

#head()

hsap.genes <- readLines('../data/hsap.genes')

sum(hsap.map.use$symbol %in% hsap.genes)

mtre <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
wtre <- readRDS('~/Documents/erc/data/saves/worms17.trees.rds')

hsap.cele.ortho <- fread('../data/inparanoid.ortholog.maps/cele_hsap/cleanedparsedOutput.C.elegans-H.sapiens') %>%
     dplyr::select(c(2,3)) %>% setnames(c('celeuniprot','hsapuniprot'))
cele.ogwbgeneuniprot <- fread('../data/worm.ogwbgeneiduniprot.map', header = T)
head(hsap.map.use)

hsap.cele.map <- left_join(cele.ogwbgeneuniprot, hsap.cele.ortho, by = 'celeuniprot') %>% na.omit() %>%
     left_join(hsap.map.use, by = c('hsapuniprot'='uniprot')) %>% na.omit() %>%
     setnames(c('og','txid','wbgnid','celeuniprot','hsapuniprot','symbol'))

length(hsap.cele.map$symbol)
length(unique(hsap.cele.map$symbol))
head(sort(table(hsap.cele.map$symbol), decreasing = T))
length(intersect(hsap.cele.map$symbol, names(mtre$trees)))

hsap.cele.map.unqhsapsymbol <- hsap.cele.map %>%
     group_by(symbol) %>%
     filter(row_number() <= 1) %>% 
     ungroup()

write.table(hsap.cele.map.unqhsapsymbol, file = '../data/cele.OgwbgnUniprotSymbol.map', 
            quote = F, sep = '\t', row.names = F)
hsap.cele.map.unqhsapsymbol <- fread('../data/cele.OgwbgnUniprotSymbol.map')
which(!hsap.cele.map.unqhsapsymbol$og %in% rownames(erc)) #=1673
hsap.cele.map.unqhsapsymbol[1673,]
write.table(hsap.cele.map.unqhsapsymbol[-1673,], file = '../data/cele.OgwbgnUniprotSymbol.map', 
            quote = F, sep = '\t', row.names = F)
write(hsap.cele.map.unqhsapsymbol$symbol[-1673], file = '../data/hsap.cele.genes')

