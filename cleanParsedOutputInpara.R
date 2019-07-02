require(dplyr)
require(data.table)
fil <- '../data/inparanoid.ortholog.maps/cele_hsap/parsedOutput.C.elegans-H.sapiens'
#fil <- '../data/inparanoid.ortholog.maps/dmel_hsap/parsedOutput.D.melanogaster-H.sapiens'
#fil <- '../data/inparanoid.ortholog.maps/hsap_scer/parsedOutput.H.sapiens-S.cerevisiae'
#fil <- '../data/inparanoid.ortholog.maps/dmel_scer/parsedOutput.D.melanogaster-S.cerevisiae'
sampdat <- fread(fil,header=F)
filn <- strsplit(fil,'/')[[1]]
filnd <- paste(strsplit(fil,'/')[[1]][1:(length(filn)-1)], collapse = '/')
filnf <- filn[length(filn)]
#head(arrange(sampdat, V6))
#filter(sampdat, V1 == 'G5EC46')
nrow(sampdat)
length(unique(sampdat$V1))
length(unique(sampdat$V2))

n2ewdat <- as.data.frame(sampdat) %>%
     filter(V7 == 100, V8 == 100) %>%
     group_by(V2) %>% arrange(desc(V7)) %>% filter(row_number() <= 1) %>%
     group_by(V3) %>% arrange(desc(V7)) %>% filter(row_number() <= 1)

n1ewdat <- as.data.frame(sampdat) %>%
     group_by(V2) %>% arrange(desc(V7)) %>% filter(row_number() <= 1)

newdat <- as.data.frame(sampdat) %>%
     group_by(V2) %>% arrange(desc(V7)) %>% filter(row_number() <= 1) %>%
     group_by(V3) %>% arrange(desc(V7)) %>% filter(row_number() <= 1)

#quantile(newdat$V6)
#quantile(newdat$V7)
#baddat <- filter(newdat, V6 < 90 | V7 < 90)
#findat <- filter(newdat, V6 >= 90 & V7 >= 90)
write.table(n2ewdat, file = paste0(filnd,'/cleaned',filnf), quote = F, col.names = F, row.names = F)
#nbad = nrow(baddat)
#}

sp1 <- strsplit(filn[4],'_')[[1]][1]
sp2 <- strsplit(filn[4],'_')[[1]][2]

write.table(n2ewdat$V2, file = paste0(filnd,'/cleaned',filnf,'.',sp1,'.uniprot.id'), quote = F, col.names = F, row.names = F)
write.table(n2ewdat$V3, file = paste0(filnd,'/cleaned',filnf,'.',sp2,'.uniprot.id'), quote = F, col.names = F, row.names = F)