args <- commandArgs(trailingOnly = TRUE)

ii <- as.numeric(args[1])

for(jj in c(1:5)){
     assign(paste0('fd',jj), readRDS(paste0('/zfs1/nclark/rap119/coevolution/data/integration/int',jj,'/int',ii,'_',jj,'.rds')))
}
fd <- cbind(fd1,fd2[,6:10],fd3[,6:10],fd4[,6:10],fd5[,6:10]) 
fd$sumnlogpv1 <- rowSums(fd[,c('nlogpv11','nlogpv12','nlogpv13','nlogpv14','nlogpv15')], na.rm = T)
fd$sumnlogpvbest <- rowSums(fd[,c('nlogpvbest1','nlogpvbest2','nlogpvbest3','nlogpvbest4','nlogpvbest5')], na.rm = T)
fd <- select(fd,c(1:5,31,32,6:30)) %>%
     arrange(desc(sumnlogpvbest))
saveRDS(fd, file = paste0('/zfs1/nclark/rap119/coevolution/data/integration/intall/int',ii,'.rds'))