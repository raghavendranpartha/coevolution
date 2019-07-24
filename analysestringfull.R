require(dplyr)
require(data.table)
require(tidyr)
require(ggplot2)

string.hg.alias.full <- fread('../data/public/stringhuman9606.protein.aliases.v11.0.txt', header = F)
string.hg.alias <-string.hg.alias.full %>% 
     mutate(V1 = gsub('9606\\.','',V1)) %>%
     filter(grepl('Ensembl_gene',V3)) %>%
     dplyr::select(c(1,2))

hg.ens.sym <- fread('../data/humanSymbolEnsembl.raw.map', header = T) %>%
     dplyr::select(c(2,4)) %>% unique()
string.hg.alias.sym <- left_join(string.hg.alias, hg.ens.sym, by = c('V2'='hg19.ensGene.name2')) %>%
     na.omit() %>% setnames(c('V1','V2','V3')) %>%
     filter(V3 %in% hsap.genes) 
rm('string.hg.alias')
rm('hg.ens.sym')
rm('string.hg.alias.full')

string.full <- fread('../data/public/stringdb/9606.protein.links.full.v11.0.newcombinedscore.txt', header = F) %>%
     mutate(protein1 = gsub('9606\\.','',V1)) %>%
     mutate(protein2 = gsub('9606\\.','',V2))
setnames(string.full,c('V1','V2','score','protein1','protein2'))
setDT(string.full)
setDT(string.hg.alias.sym)
string.full[string.hg.alias.sym,on=.(protein1=V1),g1:=i.V3]
string.full[string.hg.alias.sym,on=.(protein2=V1),g2:=i.V3]
string.full <- string.full[!is.na(g1) & !is.na(g2)]

string.full$interc <- calcpairid(string.full$g1,string.full$g2,intfterc)
string.full$merc <- calcpairid(string.full$g1, string.full$g2, merc.sym)
string.full$pair <- calcpairid(string.full$g1, string.full$g2, pairmap)
string.full$nds <- calcpairid(string.full$g1, string.full$g2, ndatafterc)

string.full.use <- filter(string.full, !g1 %in% orgzfg, !g2 %in% orgzfg,!pair %in% badblastpairs$pair) %>%
     select(c('pair','interc','merc','score','nds')) %>% filter(pair!=0,nds!=0)%>%
     group_by(pair) %>% arrange(desc(score)) %>% filter(row_number() <= 1)
write.table(string.full.use, file = '../data/public/stringdb/string.full.use.cleaned.tsv',
            quote = F, sep = '\t', row.names = F)
     

string.full.use <- fread('../data/public/stringdb/string.full.use.cleaned.tsv', header = T)


bgdbig <- readRDS('../data/toperc/bgdbig.rds')
set.seed(2)
qs.any <- quantile(string.full.use$score, c(0.98,0.99,0.995,0.998,1))
#qs.any <- c(900,950,975,990,999)
qs.any <- quantile(arrange(string.full.use,desc(score))$score[1:1000000],c(0.9,0.95,0.98,0.99,1))
string.full.controlpairs.with10pc <- ndataset.matched.controlpairs(filter(string.full.use,score > qs.any[1])$nds, nPerLevelmult = 5)

pdf('../figures/stringboxplots/string.full.interc.boxplot.pdf', height = 6, width = 8,family='FreeSans')
plotbox2(string.full.controlpairs.with10pc$interc,xx = string.full.use$score, 
         yy = string.full.use$interc, breaksi = qs.any, leglabsi = c('Control','Top 10-5%','Top 5-2%','Top 2-1%','Top 1%'), titl = paste0('STRING PPIs'), xlab = 'Confidence score', ylab = 'Integrated ERC', ylimi=c(-4,8))
dev.off()


pdf('../figures/stringboxplots/string.full.merc.boxplot.pdf', height = 6, width = 8,family='FreeSans')
plotbox2(string.full.controlpairs.with10pc$merc,xx = string.full.use$score, 
         yy = string.full.use$merc, breaksi = qs.any, leglabsi = c('Control','Top 10-5%','Top 5-2%','Top 2-1%','Top 1%'), titl = paste0('STRING PPIs'), xlab = 'Confidence score', ylab = 'Mammal ERC', ylimi=c(-4,8))
dev.off()
