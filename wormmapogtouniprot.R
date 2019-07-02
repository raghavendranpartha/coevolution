require(dplyr)
require(data.table)

wtre <- readRDS('~/Documents/erc/data/saves/worms17.trees.rds')
length(names(wtre$trees))

ogtotxidraw <- fread('../data/worms.Orthogroups.csv.1to1_FILTERED_MinSp6.txt', header = T) %>% dplyr::select(c(1,10)) %>% setnames(c('og','txid')) %>% filter(txid != '---')
#Run fol in bash
#zgrep ">" caenorhabditis_elegans.PRJNA13758.WBPS12.CDS_transcripts.fa.gz > worm.transcriptotwbgeneid.map

txidtowbgeneid <- fread('../data/worm.transcriptotwbgeneid.map', header = F) %>%
     rowwise() %>%
     mutate(V1 = gsub('>','',V1)) %>%
     mutate(V2 = gsub('gene=','',V2)) 

ogtowbgeneid <- left_join(ogtotxidraw,txidtowbgeneid, by = c('txid'='V1'))
write.table(ogtowbgeneid, file = '../data/cele.OgWbgnTxID.map', quote = F, sep = '\t', row.names = F)

# trim <- function (x) gsub("^\\s+|\\s+$", "", x)
# cenumbertouniprot <- read.fwf('../data/celegans_uniprotcenumbermap.txt', widths = c(38,7,13,26)) %>%
#      mutate(V3 = trim(V3)) %>%
#      mutate(V2 = trim(V2))
# # sum(hsap.cele.ortho$celeuniprot %in% cenumbertouniprot$V3)

# hsap.cele.ortho <- fread('../data/inparanoid.ortholog.maps/cele_hsap/cleanedparsedOutput.C.elegans-H.sapiens') %>%
#      select(c(2,3)) %>% setnames(c('celeuniprot','hsapuniprot'))
cele.wbg.uniprot <- fread('../data/celegans.wbgeneiduniprot.parasitemart.txt', header = T) %>% select(c(4,5)) %>% unique() %>%
     setnames(c('celeuniprot','wbgeneid')) %>% filter(celeuniprot != '')
length(unique(cele.wbg.uniprot$wbgeneid))
ogtowbgeneiduniprot <- left_join(ogtowbgeneid,cele.wbg.uniprot, by = c('V2'='wbgeneid')) %>%
     na.omit()
length(unique(ogtowbgeneiduniprot$V2))



write.table(ogtowbgeneiduniprot, file = '../data/worm.ogwbgeneiduniprot.map', quote = F, sep = '\t',
            row.names = F)
# sum(hsap.cele.ortho$celeuniprot %in% ogtowbgeneiduniprot$celeuniprot)
# ogtowbgeneiduniprot.use <- left_join(ogtowbgeneiduniprot, hsap.cele.ortho, by = 'celeuniprot') %>%
#      na.omit()
