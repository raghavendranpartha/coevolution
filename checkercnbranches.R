library(RERconverge)
source('~/Documents/RERc/RERconverge/R/RERfuncs.R')
source('~/Documents/RERc/RERconverge/R/RcppExports.R')
source('~/Documents/RERc/RERconverge/R/projection_coevo.R')
source('~/Documents/erc/code/funcsERCfromProjections.R')

gns <- c('TSN','MYL12A')
mtre <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.rds')
mwts <- readRDS('~/Documents/erc/data/saves/mamm63nt.trees.weights.rds')
mrers <- readRDS('~/Documents/rermethods/data/mamm63nt.trees.scaledrers.sqrt.wt.rds')

correlateERCGeneListTrees=function(treesObj, genelist, cutoff=NULL, transform="sqrt", weighted=T, weights = NULL, useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05){
     
          totalpairs <- choose(treesObj$numTrees,2)
          treesObj = mtre
          genelist <- gns
          cutoff<-NULL
          transform="sqrt"
          weighted=T
          weights = mwts
          useSpecies=NULL
          min.sp=5
          scale=T
          doOnly=NULL
          maxT=NULL
          scaleForPproj=F
          mean.trim=0.05
     
     totalpairs <- choose(length(genelist),2)
     
     if(is.null(cutoff)){
          cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
          message(paste("cutoff is set to", cutoff))
     }
     if (weighted){
          if(is.null(weights)){
               message("Weights can not be NULL if weighted = T; computing weights")
               weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
          }
          residfunc=fastLmResidMatWeighted
     }else{
          residfunc=fastLmResidMat
     }
     
     if (is.null(useSpecies)){
          useSpecies=treesObj$masterTree$tip.label
          #mappedEdges=trees$mappedEdges
     }
     if(is.null(maxT)){
          maxT=treesObj$numTrees
     }
     if(transform!="none"){
          transform=match.arg(transform,c("sqrt", "log"))
          transform=get(transform)
     }else{
          transform=NULL
     }
     #cm is the names of species that are included in useSpecies and the master tree
     cm=intersect(treesObj$masterTree$tip.label, useSpecies)
     sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
     if (length(sp.miss) > 0) {
          message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                       collapse = ",")))
          
     }
     
     corout=matrix(NA,treesObj$numTrees,treesObj$numTrees)
     diag(corout) = 1
     
     geneis <- sort(match(genelist, names(treesObj$trees)))
     #genejs <- match(genelist, names(treesObj$trees))
     
     maxn2=treesObj$report%*%t(treesObj$report)
     maxn2[upper.tri(maxn2, diag = T)] = NA
     corout[which(maxn2 == 0)] <- Inf
     
     
     #      if(is.null(doOnly)){
     #           doOnly=1
     #      }else{
     #           maxT=1
     #      }
     pairsdone <- 0
     print(maxT)
     for (i in geneis){
          for(j in geneis){
               if(i < j){
                    next
               }
               print(paste0(i,'',j))
               #print(j)
               #i <- 480
               #j <- 2
               i <- 7265
               j <- 1715
               if(is.na(corout[i,j])){
                    #tree1=treesObj$trees[[i]]
                    #get the common species, prune and unroot
                    both=intersect(intersect(treesObj$trees[[i]]$tip, treesObj$trees[[j]]$tip), cm)
                    nb=length(both)
                    if(nb>2){
                         tree1=unroot(pruneTree(treesObj$master,both))
                    }
                    #find all the genes that contain all of the species in tree1
                    allreport=treesObj$report[,both]
                    if(nb>1){
                         ss=rowSums(allreport)
                    }else{
                         ss=allreport
                    }
                    iiboth=which(ss==length(both))
                    ai2=which(maxn2[iiboth, iiboth]==nb, arr.ind = T)
                    
                    if(length(both)<min.sp){
                         corout[matrix(iiboth[ai2],nrow=nrow(ai2))]=Inf
                         ai2rev <- ai2
                         ai2rev[,1] <- ai2[,2]
                         ai2rev[,2] <- ai2[,1]
                         corout[matrix(iiboth[ai2rev],nrow=nrow(ai2rev))]=NA
                         pairsdone <- pairsdone+nrow(ai2)
                         #message(paste0(pairsdone," done out of ",totalpairs))
                         next
                    }    
                    
                    
                    if(T){
                         
                         
                         ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
                         
                         ii= treesObj$matIndex[ee[, c(2,1)]]
                         
                         allbranch=treesObj$paths[iiboth,ii]
                         if(weighted){
                              allbranchw=weights[iiboth,ii]
                         }
                         if(scaleForPproj){
                              nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
                         }else{
                              nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
                         }
                         
                         iibad=which(allbranch<cutoff)
                         #don't scale
                         #allbranch=scaleMat(allbranch)
                         if(!is.null(transform)){
                              nv=transform(nv)
                              allbranch=transform(allbranch)
                         }
                         allbranch[iibad]=NA
                         
                         if(!scale){
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                                   
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv), allbranchw)
                                   
                              }
                         }else{
                              
                              if(!weighted){
                                   proj=residfunc(allbranch, model.matrix(~1+nv))
                              }else{
                                   
                                   proj=residfunc(allbranch, model.matrix(~1+nv),allbranchw)
                              }
                              
                              proj=scale(proj, center = F)
                              
                         }
                         
                         
                         #we have the projection
                         
                         ai2[which(ai2[,1]==1204),]
                         
                         for(m in 1:nrow(ai2)){
                              m = 293
                              #print(m)
                              k=sort(ai2[m,])[1]
                              l=sort(ai2[m,])[2]
                              #tmpcor=cor(proj[k,], proj[l,], method = 'p', use = 'pair')
                              re1 <- proj[k,]
                              re2 <- proj[l,]
                              print(re1[!is.na(re1)])
                              print(re2[!is.na(re2)])
                              goodinds <- !is.na(re1)&!is.na(re2)
                              if(sum(length(goodinds)) >= min.sp){
                                   plot(re1[goodinds],re2[goodinds])
                              }
                              if(sum(!is.na(proj[k,])&!is.na(proj[l,])) > 4){
                                   cres=cor.test(win(proj[k,],3), win(proj[l,],3), method='pearson', exact=F)
                                   corout[iiboth[l], iiboth[k]]=cres$est
                              }else{
                                   corout[iiboth[l], iiboth[k]]=Inf
                              }
                              corout[iiboth[k], iiboth[l]]=1
                         }
                         pairsdone <- pairsdone+nrow(ai2)
                         #message(paste0(pairsdone," done out of ",totalpairs))
                    }
                    
               }
               
          }
     }
     rownames(corout) <- names(treesObj$trees)
     colnames(corout) <- names(treesObj$trees)
     corout[geneis,geneis]
}