args <- commandArgs(trailingOnly = T)

#infol <- '../data/geneSimBP/MF/onembgdpairs.bestbitscorepclt0p1.tsv'
infol <- args[1]
onttype = args[2]
outf <- tail(strsplit(infol,'/')[[1]],1)
resfol <- paste0('../data/geneSim/',onttype,'/',outf)
#print(getwd())
#print(resfol)
#print(dir.exists(resfol))
fs <- dir(resfol)
#print(fs)
fs.vecs <- lapply(fs, function(x) readRDS(file.path(resfol,x)))
vec<-unlist(fs.vecs)
#print(vec[1:10])

saveRDS(vec, paste0('../data/geneSim/',onttype,'/',outf,'.rds'))
#v2 <- readRDS('../data/geneSimBP/MF/onembgdpairs.bestbitscorepclt0p1.tsv.rds')
#identical(v2,vec)