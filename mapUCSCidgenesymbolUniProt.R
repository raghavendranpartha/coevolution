require(dplyr)
require(data.table)
require(tidyr)
library(RERconverge)

dmel.uniprot.whsaportholog <- readLines('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens.dmel.uniprot.id')
scer.uniprot.whsaportholog <- readLines('../data/inparanoid.ortholog.maps/hsap_scer/cleanedparsedOutput.H.sapiens-S.cerevisiae.scer.uniprot.id')

hg19.uniprot.wdmelortholog <- readLines('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens.hsap.uniprot.id')
hg19.uniprot.wscerortholog <- readLines('../data/inparanoid.ortholog.maps/hsap_scer/cleanedparsedOutput.H.sapiens-S.cerevisiae.hsap.uniprot.id')

hg19.uniprot.symbol.all <- fread('../data/genesymbolUCSCnameUniProt.raw.map', header = T)
head(hg19.uniprot.symbol.all)

with(hg19.uniprot.symbol.all, identical(hg19.knownGene.proteinID, hg19.kgXref.spID))
nonmatch <- with(hg19.uniprot.symbol.all, which(hg19.knownGene.proteinID != hg19.kgXref.spID))

head(hg19.uniprot.symbol.all[nonmatch,])

hg19.uniprot.symbol.dmelo.use <- hg19.uniprot.symbol.all %>% filter(hg19.kgXref.spID %in% hg19.uniprot.wdmelortholog) %>%
     dplyr::select(c(3,5)) %>% unique()

hg19.uniprot.symbol.scero.use <- hg19.uniprot.symbol.all %>% filter(hg19.kgXref.spID %in% hg19.uniprot.wscerortholog) %>%
     dplyr::select(c(3,5)) %>% unique()

length(unique(hg19.uniprot.symbol.dmelo.use$hg19.kgXref.geneSymbol))
length(unique(hg19.uniprot.symbol.scero.use$hg19.kgXref.geneSymbol))

dtres <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')
#ytres <- readTrees('~/Documents/yeastERC/data/all_Scer.tre')
#saveRDS(ytres, file = '~/Documents/erc/data/saves/all_Scer.tre.rds')
ytres <- readRDS('~/Documents/erc/data/saves/all_Scer.tre.rds')

head(names(ytres$trees))

# v2hg19.uniprot.symbol.dmelo.use <- hg19.uniprot.symbol.all %>% filter(hg19.knownGene.proteinID %in% hg19.uniprot.wdmelortholog) %>%
#      dplyr::select(c(2,5)) %>% unique()
# 
# v2hg19.uniprot.symbol.scero.use <- hg19.uniprot.symbol.all %>% filter(hg19.knownGene.proteinID %in% hg19.uniprot.wscerortholog) %>%
#      dplyr::select(c(2,5)) %>% unique()
# 
# length(v2hg19.uniprot.symbol.scero.use$hg19.kgXref.geneSymbol %in% hg19.uniprot.symbol.scero.use$hg19.kgXref.geneSymbol)
# length(v2hg19.uniprot.symbol.dmelo.use$hg19.kgXref.geneSymbol %in% hg19.uniprot.symbol.dmelo.use$hg19.kgXref.geneSymbol)


#droso map

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
dros.fbg.uniprot <- read.fwf('../data/fly_uniprotFBGN.map.txt', widths = c(40,6,6,11), header = F, stringsAsFactors = F) %>%
     dplyr::select(c(2,4)) %>%
     mutate(V2 = trim(V2))

scer.fbg.uniprot <- read.fwf('../data/yeast_uniprot.map.txt', widths = c(75,9,11,6,40), header = F, stringsAsFactors = F) %>%
     dplyr::select(c(2,4)) %>%
     mutate(V2 = trim(V2)) %>%
     filter(V2 %in% names(ytres$trees))

sum(scer.fbg.uniprot$V4 %in% scer.uniprot.whsaportholog)
#sum(scer.fbg.uniprot$V2 %in% names(ytres$trees))


sum(dmel.uniprot.whsaportholog %in% dros.fbg.uniprot$V2)
sum(scer.uniprot.whsaportholog %in% scer.fbg.uniprot$V4)

sum(filter(dros.fbg.uniprot, V2 %in% dmel.uniprot.whsaportholog)$V2 %in% names(dtres$trees))
sum(filter(scer.fbg.uniprot, V4 %in% scer.uniprot.whsaportholog)$V2 %in% names(ytres$trees))


org.Dm.egUNIPROT



