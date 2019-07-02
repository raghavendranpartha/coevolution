require(dplyr)
require(data.table)

dtres <- readRDS('~/Documents/erc/data/saves/alldroso22.tre.rds')
# 
# dmelidmap <- fread('../data/dmelorthofinder_1to1_FILTERED.txt', header = F) %>%
#      select(c(V1,V15)) %>% rowwise() %>%
#      mutate(fbgn = strsplit(V15,'_')[[1]][1])
# 
# dmelidmapgood <- filter(dmelidmap, fbgn != '---')
dmelorig <- fread('../data/dmelortho2gene.txt', header = F)

#identical(dmelidmapgood$fbgn, dmelorig$V2)

dtres$numTrees
length(unique(names(dtres$trees)))

head(sort(names(dtres$trees)))
tail(sort(names(dtres$trees)))

fbgunip <- fread('../data/FlyBase_Fields_download.txt', skip = 4, header = F, fill = T) %>% select(c(3,6)) %>%
     unique()

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
fbgunipv2 <- read.fwf('../data/fly_uniprotFBGN.map.txt', widths = c(40,6,6,11), header = F, stringsAsFactors = F) %>%
     dplyr::select(c(4,2)) %>%
     mutate(V2 = trim(V2)) %>%
     mutate(V4 = trim(V4)) %>%
     filter(V4 != '')

sum(fbgunipv2$V4 %in% dmelorig$V2)

length(unique(fbgunip.tomap$V3))

dmelwhsapunip <- fread('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens.dmel.uniprot.id', header = F)
dmelwscerunip <- fread('../data/inparanoid.ortholog.maps/dmel_scer/cleanedparsedOutput.D.melanogaster-S.cerevisiae.dmel.uniprot.id', header = F)
alldmel <- unique(c(dmelwscerunip$V1,dmelwhsapunip$V1))
#dmelwhsapunip.st <- sapply(strsplit(dmelwhsapunip$V1,''),'[[',1)
#dmelwhsapunip.len <- sapply(strsplit(dmelwhsapunip$V1,''), length)


fbgunipv2.tomap <- filter(fbgunipv2, V2 %in% alldmel) %>%
     setnames(c('fbgn','uniprot'))
fbgunip.tomap <- fbgunip %>% filter(V6 %in% alldmel) %>%
     setnames(c('fbgn','uniprot'))

length(intersect(fbgunip.tomap$fbgn, fbgunipv2.tomap$fbgn))

checktwomaps <- left_join(fbgunipv2.tomap, fbgunip.tomap, by = 'fbgn') %>%
     na.omit()
with(checktwomaps, identical(uniprot.x,uniprot.y))

combined.fbgunip.map <- rbind(fbgunip.tomap,
                              filter(fbgunipv2.tomap, !fbgn %in% fbgunip.tomap$fbgn)) %>%
     left_join(dmelorig, by = c('fbgn'='V2')) %>% na.omit() %>%
     setnames(c('fbgn','uniprot','og','symbol'))

write.table(combined.fbgunip.map, file = '../data/dmel.OgFbgnUniprotSymbol.map',
            quote = F, sep = '\t', row.names = F, col.names = T)
