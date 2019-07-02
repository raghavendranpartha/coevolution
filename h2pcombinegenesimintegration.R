#fn <- '../data/toperc/int.sumnlogpv1gteq2paircgteq3bitscorepclt10.tsv'
fn <- '../data/toperc/int.sumnlogpv1lt2paircgteq3bitscorepclt10.tsv'
int.min2.pairc3.dup <- fread(fn, header = T)

outf <- gsub('tsv','withgeneSim.tsv',fn) 
     
alldf <- lapply(1:1000, function(x){
     print(x)
     readRDS(paste0('../data/geneSimBP/int.sumnlogpv1lt2paircgteq3bitscorepclt10.',x,'.rds'))
})
int.min2.pairc3.dup$geneSimbp <- do.call(c, alldf)
int.min2.pairc3.dup <- dplyr::select(int.min2.pairc3.dup, c(1:4,35,5:34))
write.table(int.min2.pairc3.dup, file = outf,
            quote = F, sep = '\t', row.names = F)