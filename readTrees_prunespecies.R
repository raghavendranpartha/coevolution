readTrees_prunespp=function(file, spp = NULL, minsp = 5, max.read=NA){
#       file = '~/Documents/Nematodes/worms25.trees'
#       spp <- c('DipCor','DipPac','TelCir','OesDen','AncDuo','AncCan','CaeAng')
#       max.read=NA
#       minsp = 5
     tmp=scan(file, sep="\t", what="character")
     trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
     treenames=character()
     maxsp=0; # maximum number of species
     
     for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
          if (i %% 2==1){
               treenames=c(treenames, tmp[i])
          }
          else{
               #print(i)
               if(is.null(spp)){
                    tmpt <- read.tree(text=tmp[i])
                    if(length(tmpt$tip) >= minsp){
                         trees[[i/2]]=unroot(tmpt)
                    }else{
                         trees[[i/2]] <- NA
                    }
               }else{
                    tmpt <- drop.tip(read.tree(text=tmp[i]), spp)
                    if(length(tmpt$tip) >= minsp){
                         trees[[i/2]]=unroot(tmpt)
                    }else{
                         trees[[i/2]] <- NA
                    }
               }
               
               #check if it has more species
               if(!is.na(trees[[i/2]])){
                    if(length(trees[[i/2]]$tip.label)>maxsp){
                         maxsp=length(trees[[i/2]]$tip.label)
                         allnames=trees[[i/2]]$tip.label
                    }
               }
          }
     }
#      nsp <- sapply(1:length(trees), function(x){
#           length(trees[[x]]$tip)
#      })
     validtrees <- which(!is.na(trees))
     trees <- trees[validtrees]
     names(trees)=treenames[validtrees]
     treesObj=vector(mode = "list")
     treesObj$trees=trees
     treesObj$numTrees=length(trees)
     treesObj$maxSp=maxsp
     
     message(paste("max is ", maxsp))
     
     report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
     colnames(report)=allnames
     rownames(report)=treenames[validtrees]
     print('here')
     for ( i in 1:nrow(report)){
          ii=match(allnames, trees[[i]]$tip.label)
          report[i,]=1-is.na(ii)
          
     }
     treesObj$report=report
     
     
     
     ii=which(rowSums(report)==maxsp)
     
     #Create a master tree with no edge lengths
     master=trees[[ii[1]]]
     master$edge.length[]=1
     treesObj$masterTree=master
     
     
     
     
     treesObj$masterTree=rotateConstr(treesObj$masterTree, sort(treesObj$masterTree$tip.label))
     #this gets the abolute alphabetically constrained order when all branches
     #are present
     tiporder=treeTraverse(treesObj$masterTree)
     
     #treesObj$masterTree=CanonicalForm(treesObj$masterTree)
     
     for ( i in 1:treesObj$numTrees){
          #print(i)
          treesObj$trees[[i]]=rotateConstr(treesObj$trees[[i]], tiporder)
          
     }
     
     
     
     ap=allPaths(master)
     treesObj$ap=ap
     matAnc=(ap$matIndex>0)+1-1
     matAnc[is.na(matAnc)]=0
     
     paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
     for( i in 1:treesObj$numTrees){
          print(i)
          paths[i,]=allPathMasterRelative(treesObj$trees[[i]], master, ap)
     }
     paths=paths+min(paths[paths>0], na.rm=T)
     treesObj$paths=paths
     treesObj$matAnc=matAnc
     treesObj$matIndex=ap$matIndex
     treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))
     
     ii=which(rowSums(report)==maxsp)
     if(length(ii)>20){
          message (paste0("estimating master tree branch lengths from ", length(ii), " genes"))
          tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})
          
          allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
          allEdge=scaleMat(allEdge)
          allEdgeM=apply(allEdge,2,mean)
          treesObj$masterTree$edge.length=allEdgeM
     }
     else{
          message("Not enough genes with all species present: master tree has no edge.lengths")
     }
     colnames(treesObj$paths)=namePathsWSpecies(treesObj$masterTree)
     class(treesObj)=append(class(treesObj), "treesObj")
     treesObj
}