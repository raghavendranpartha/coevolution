#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

#library(org.Hs.eg.db)

require(dplyr)
require(data.table)

hsap.genes <- readLines('../data/hsap.genes')

#genemap <- fread('../data/Supplementary Table 3. Predicted PPI across 122 species.txt')
#densids <- unlist(strsplit(genemap$HsPPIEnsembl, '---'))
#dgenes <- unlist(strsplit(genemap$HsPPIGene, '---'))
#ensgenemap <- data.frame(ensid = densids, sym = dgenes, stringsAsFactors = F)

geneomaage <- fread('../data/OMAproteinage.gene')
# biom <- fread('../data/mart_export.txt') %>%
#      setnames(c('ensgene','enstrans','hgncsym','hgncid'))
ens.sym <- fread('../data/humanSymbolEnsembl.raw.map') %>%
     setnames(c('enstid1','ensgene','enstid2','symbol'))
#ensucsc <- fread('../data/humanucscIdEnsembl.raw.map', header = T) %>%
#     filter(hg19.knownToEnsembl.name != 'n/a')

# sum(hsap.genes %in% biom$hgncsym)
# sum(geneomaage$Gene %in% filter(biom, hgncsym %in% hsap.genes)$ensgene)

ens.sym.use <- filter(ens.sym, symbol %in% hsap.genes) %>%
     select(c(2,4)) %>% unique()
ens.sym.use.vec <- ens.sym.use$symbol
names(ens.sym.use.vec) <- ens.sym.use$ensgene

geneomaage <- geneomaage %>%
     mutate(symbol = ens.sym.use.vec[Gene]) %>% na.omit()

#read in cplx.hsap from cleanmarcottecomplexes.R
marc.cplx.hsap <- unique(unlist(cplx.hsap, use.names = F))     

cplx.age <- sapply(1:length(cplx.hsap), function(x){
     mean(filter(geneomaage, symbol %in% cplx.hsap[[x]])$ProteinAge == 'new')
})
quantile(cplx.age, probs = seq(0,1,length.out = 11))

write.table(geneomaage, file = '../data/Marcotte_981complexes.genes.omaage.tsv',
            quote = F, row.names = F, sep = '\t')

# sum(hsap.genes %in% ens.sym$symbol)
# sum(geneomaage$Gene %in% filter(ens.sym, symbol %in% hsap.genes)$ensgene)



#sum(geneomaage$Gene %in% biom$ensgene)
#sum(geneomaage$Gene %in% ens.sym$ensgene)
#sum(geneomaage$Gene %in% ensucsc$hg19.ensGene.name2)

#hsap.genes.ucid <- fread('~/Documents/mammal_data/data/ucidCorrect_commonName.map')
#sapply(hsap.genes.ucid$ucid, function(x) sum(grepl(x, ensucsc$hg19.knownToEnsembl.name)))

#sum(geneomaage$Gene %in% ensucsc$hg19.ensGene.name2)
sum(geneomaage$Gene %in% hsapenstosym$ensid)

# mts <- sapply(hsap.genes.ucid$ucid, function(x){
#      length(grep(x, ensucsc$hg19.knownToEnsembl.name, value = T))
# })
# sum(mts == 1)
# sum(mts == 0)
# table(mts)






# hsapenstosym <- fread('../data/humanSymbolEnsembl.raw.map', header = T) %>%
#      select(c(2,4)) %>% unique() %>%
#      filter(hg19.ensemblToGeneName.value %in% hsap.genes) %>%
#      setnames(c('ensid','symbol'))

#length(unique(hsapenstosym$symbol))

#sum(hsap.genes %in% hsapenstosym$hg19.ensemblToGeneName.value)

#sum(hsap.genes %in% hsapenstosym$hg19.ensemblToGeneName.value)
#noma <- setdiff(hsap.genes,hsapenstosym$hg19.ensemblToGeneName.value)

cplx <- strsplit(dat$GeneName,';')
allens <- unlist(strsplit(dat$EnsemblID,';'))
#datallgenes <- unlist(cplx)
cplx.hsap <- sapply(cplx, function(x){
     intersect(x,hsap.genes)
})
datallgenes <- unique(unlist(cplx.hsap))
sum(datallgenes %in% hsap.genes)
sum(datallgenes %in% hsapenstosym$symbol)


length(unique(hsapenstosym$hg19.ensemblToGeneName.value))
length(unique(hsapenstosym$hg19.ensGene.name2))

sum(table(hsapenstosym$hg19.ensemblToGeneName.value) > 1)
sum(geneomaage$Gene %in% hsapenstosym$ensid)
length(unique(datallgenes))
#sum(unique(datallgenes) %in% hsap.genes)
#sum(unique(datallgenes) %in% hsapenstosym$symbol)
sum(unique(allens) %in% hsapenstosym$ensid)

