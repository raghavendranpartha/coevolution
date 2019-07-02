require(dplyr)
require(data.table)
require(tidyr)

#get hg19 genes from addspalaxmamm62/code
#makeblastdb -in allgeneshg19.fa -parse_seqids -dbtype nucl -out fadb/allgeneshg19
#blastn -query allgeneshg19.fa -db fadb/allgeneshg19 -outfmt 6 -perc_identity 50 > allgeneshg19blast.fmt6.out
#blastn -query allgeneshg19.fa -db fadb/allgeneshg19 -word_size 11 -outfmt 5 -perc_identity 50 > allgeneshg19blast.fmt5.raw.ws7pid50.out
#perl blast_parser.pl 10 allgeneshg19blast.fmt5.raw.ws7pid50.out > allgeneshg19blast.fmt5.ws7pid50.preproc.out
#cut -f1-9 allgeneshg19blast.fmt5.ws7pid50.preproc.out > allgeneshg19blast.fmt5.ws7pid50.out

blast6 <- fread('~/Documents/addspalaxmamm62/data/allgeneshg19blast.fmt6.raw.ws7pid50.out', header = F)
setnames(blast6,c('gene1','gene2','perc_id','aln.len','numMismatch','numGapOpen','aln.start.query','aln.end.query','aln.start.subj','aln.end.subj','evalue','bitscore'))

View(filter(blast6, V12 < 50))

genepair <-'ZFP90|ZNF189'
genepair <-'A1BG|IGSF1'
fd <- filter(blast6, grepl(genepair, V1),grepl(genepair, V1)) %>%
     rowwise() %>% filter(V1 != V2) %>%
     mutate(pair = paste(sort(c(V1,V2)), collapse = ',')) %>% dplyr::select(c(13,1:12))
setnames(fd,c('pair','gene1','gene2','perc_id','aln.len','numMismatch','numGapOpen','aln.start.query','aln.end.query','aln.start.subj','aln.end.subj','evalue','bitscore'))
View(fd)

allpairs.blast.bitscorepc <- fread('../data/allgeneshg19.pairbitscore.tsv', header = T)
filter(allpairs.blast.bitscorepc, grepl('ZFP90,ZNF189', pair))
View(filter(allpairs.blast.bitscorepc, grepl('A1BG,IGSF1', pair)))
View(filter(allpairs.blast.bitscorepc, grepl('A2M,PZP', pair)))

#blast 5 pre-processing
blast5 <- fread('~/Documents/addspalaxmamm62/data/allgeneshg19blast.fmt5.ws7pid50.out', header = F)
fd2 <- filter(blast5, grepl(genepair, V1),grepl(genepair, V1))
fd2 <- filter(blast5, grepl('A2M|PZP', V1),grepl('A2M|PZP', V2))



gns6 <- with(blast6, c(V1, V2))
gns5 <- with(blast5, c(V1, V2))

nmatch6 <- sort(table(gns6), decreasing = T)
nmatch5 <- sort(table(gns5), decreasing = T)
singlegenes6 <- names(nmatch6)[nmatch6 == 2]
singlegenes5 <- names(nmatch5)[nmatch5 == 2]

length(intersect(singlegenes5,singlegenes6))

genesnotin5 <- setdiff(singlegenes5, singlegenes6)

multiblastout6 <- filter(blast6, ! V1 %in% singlegenes6) %>%
     rowwise() %>%
     filter(V1 != V2) %>%
     mutate(pair = paste(sort(c(V1,V2)), collapse = ',')) %>%
     separate(pair,c('g1','g2'),',',remove = F) %>%
     dplyr::select(c(13,14,15,3:12))  %>%
     
multiblastout6.final <- multiblastout6 %>% group_by(pair,bitscore) %>%
     filter(row_number() <= 1)

View(filter(blast5, grepl('A1BG|TTN', V1),grepl('A1BG|TTN', V2)))

blast5.same <- blast5 %>%
     rowwise() %>%
     filter(V1 == V2)
genebitscorewithitself <- blast5.same$V3
names(genebitscorewithitself) <- blast5.same$V2

blast5.diff <- blast5 %>% rowwise() %>%
     filter(V1 != V2) %>%
     mutate(pair = paste(sort(c(V1,V2)), collapse = ',')) %>%
     separate(pair,c('g1','g2'),',',remove = F) %>%
     dplyr::select(c(10,11,12,3:9))

blast5.diff.unq <- group_by(blast5.diff, pair) %>%
     filter(row_number() <= 1)

blast5.diff.unq.final <- blast5.diff.unq %>%
     rowwise() %>%
     mutate(bitscorepcquery = V3*1.0/genebitscorewithitself[g1]) %>%
     mutate(bitscorepcsubj = V3*1.0/genebitscorewithitself[g2]) %>%
     mutate(bestbitscorepc = max(bitscorepcsubj,bitscorepcquery))

setnames(blast5.diff.unq.final, c('pair','g1','g2','bitscore','query.length','subj.length','longestmatch.query','longestmatch.subj','totalmatch.query','totalmatch.subj', 'bitscorepc.query','bitscorepc.subj','bestbitscorepc'))
write.table(blast5.diff.unq.final, file = '../data/allgeneshg19.pairbitscore.tsv',
            quote = F, sep = '\t', row.names = F)

quantile(blast5.diff.unq.final$bestbitscorepc)
nrow(filter(blast5.diff.unq.final, bestbitscorepc >= 0.2))
badpairs <- filter(blast5.diff.unq.final, bestbitscorepc >= 0.1)
setnames(badpairs, c('pair','g1','g2','bitscore','query.length','subj.length','longestmatch.query','longestmatch.subj','totalmatch.query','totalmatch.subj',
                     'bitscorepc.query','bitscorepc.subj','bestbitscorepc'))
write.table(badpairs, file = '../data/allgeneshg19.badpairs.bitscore.tsv',
            quote = F, sep = '\t', row.names = F)

length(unique(blast5.diff$pair))
head(sort(table(blast5.diff$pair), decreasing = T))

multiblastout5 <- filter(blast5, ! V1 %in% singlegenes5) %>%
     rowwise() %>%
     mutate(pair = paste(sort(c(V1,V2)), collapse = ',')) %>%
     separate(pair,c('g1','g2'),',',remove = F) %>%
     dplyr::select(c(10,11,12,3:9)) %>%
     group_by(pair) %>%
     filter(row_number() <= 1)
multiblastout5.diff <- multiblastout5 %>%
     rowwise() %>% filter(g1 != g2)
multiblastout5.final <- multiblastout5.diff %>%
     mutate(totalmatchquerypc = totalmatchlengthquery/querylen) %>%
     mutate(totalmatchhitpc = totalmatchlengthhit/hitlen) %>%
     rowwise() %>%
     mutate(bestmatchpc = max(totalmatchquerypc, totalmatchhitpc))

length(unique(multiblastout6$pair))
pairdiff <- setdiff(unique(multiblastout6$pair), multiblastout5.diff$pair)

head(multiblastout5)
head(multiblastout6)

setnames(multiblastout5.diff, c('pair','gene1','gene2','bitscore','querylen','hitlen','longestsegmentquery','longestsegmenthit','totalmatchlengthquery','totalmatchlengthhit'))
setnames(multiblastout6,c('pair','gene1','gene2','perc_id','aln.len','numMismatch','numGapOpen','aln.start.query','aln.end.query','aln.start.subj','aln.end.subj','evalue','bitscore'))
filter(multiblastout5.final, grepl('A2M',pair))
filter(multiblastout6.final, grepl('A2M',pair))
filter(multiblastout5.final, grepl('ZXDB,ZXDC',pair))
filter(multiblastout6.final, grepl('ZXDB,ZXDC',pair))
filter(multiblastout5.final, grepl('ABCA10,ABCA8',pair))
#filter(multiblastout6, grepl('ABCA10,ABCA8',pair))
filter(multiblastout6.final, grepl('ABCA10,ABCA8',pair))


write.table(multiblastout5.final, file = '../data/allgeneshg19.blast.outfmt5.final.tsv',
            quote = F, row.names = F, sep = '\t')
write.table(multiblastout6.final, file = '../data/allgeneshg19.blast.outfmt6.final.tsv',
            quote = F, row.names = F, sep = '\t')


#filter(multiblastout6, grepl(pair, genesnotin5[1]))

multiblastout.best <- group_by(multiblastout, pair) %>%
     arrange(desc(V12)) %>%
     filter(row_number() <= 1)
     
View(filter(multiblastout.best, grepl('A2M|PZP',pair)))

multiblastout.best.diff <- multiblastout.best %>%
     rowwise() %>%
     filter(V1!=V2) %>%
     separate(pair,c('g1','g2'),',',remove = F) %>%
     dplyr::select(c(14,15,3:13))
multiblastout.best.diff$g1len <- genelen[multiblastout.best.diff$g1]
multiblastout.best.diff$g2len <- genelen[multiblastout.best.diff$g2]
multiblastout.best.diff <- multiblastout.best.diff %>%
     rowwise() %>%
     mutate(alnpc = V4/min(g1len,g2len)) %>%
     dplyr::select(c(13,1,2,14,15,16,3:12))

head(multiblastout.best.diff)

with(multiblastout.best.diff, length(unique(c(V1, V2))))

multiblastout.same <- multiblastout %>%
     rowwise() %>%
     filter(V1 == V2)

blastout <- fread('~/Documents/addspalaxmamm62/data/allgeneshg19blastv2.out',header = F) %>%
     rowwise() %>%
     filter(V1 == V2)

nmatch <- sort(table(blastout$V1), decreasing = T)
sum(nmatch == 1)
singlegenesv2 <- names(nmatch)[nmatch == 1]

multiblastoutv2 <- filter(blastout, ! V1 %in% singlegenesv2) %>%
     mutate(querymatchpc = V8*1.0/V4) %>%
     mutate(hitmatchpc = V9*1.0/V5) %>%
     rowwise() %>%
     mutate(pair = paste(sort(c(V1,V2)), collapse = '_'))
     
View(filter(multiblastout, grepl('A2M|PZP',pair)))


multiblastout.same <- multiblastout %>%
     rowwise() %>%
     filter(V1 == V2)

head(multiblastout)

with(multiblastout, sum( V4 == V8))

View(filter(multiblastout, V1 %in% c('A2M','PZP')))

multiblastout.use <- multiblastout %>%
     rowwise() %>% filter(V1 != V2) %>%
     group_by(pair) %>%
     summarise(matchpc = max(querymatchpc, hitmatchpc))

write.table(multiblastout.use, file = '../data/duplicatedgenepairs.pc.tsv',
            quote = F, row.names = F, sep = '\t')

quantile(multiblastout.use$matchpc, probs = seq(0,1,0.1))