
trn <- c('mtre','vtre','dtre','wtre','ytre')
lbls <- c('Mammal','Vertebrate','Fly','Worm','Yeast')
par(mfrow = c(1,5), mar = c(2,2,2,2))
for(ii in c(1:5)){
#     ii <- 1
     tr <- get(trn[ii])$master
     print(tr$tip)
     if(ii < 3){
          tr <- pruneTree(tr,c('Human',sample(tr$tip.label,19)))
     }
     tr$edge.length <- rep(1,length(tr$edge.length))
     plot.phylo(tr,font=2,main=lbls[ii],cex=1.25,cex.main=2,
                x.lim = c(0,16))
}
