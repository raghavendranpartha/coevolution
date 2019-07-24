
# source("https://bioconductor.org/biocLite.R")
# biocLite('STRINGdb')
require(STRINGdb)
require(dplyr)
require(data.table)


string.full.mod <- string.full
string.full.mod$cooccurence <- 0
string.full.mod$neighborhood <- 0
string.full.mod$neighborhood_transferred <- 0
string.full.mod$homology <- 0
string.full.mod$fusion <- 0
write.table(string.full.mod[,1:16], file = '../data/public/stringdb/9606.protein.links.full.v11.0.modified.txt',
            quote = F, row.names = F, sep = '\t')

stringdb <- STRINGdb$new(version = '10', species = 9606, score_threshold = 150)
intall <- stringdb$get_interactions()

string.type <- fread('../data/public/9606.protein.actions.v11.0.txt', header = T) %>% rowwise() %>%
     mutate(pair = paste(sort(c(item_id_a,item_id_b)), collapse = '_'))

string.type <- arrange(string.type, pair)

string.detailed <- fread('../data/public/9606.protein.links.detailed.v11.0.txt', header = T) %>% rowwise() %>%
     mutate(pair = paste(sort(c(protein1,protein2)), collapse = '_'))

setDT(string.detailed)
setkey(string.detailed,pair)

setDT(string.type)
setkey(string.type,pair)

head(string.type, 50)

string.full <- fread('../data/public/9606.protein.links.full.v11.0.txt', header = T) %>% rowwise() %>%
     mutate(pair = paste(sort(c(protein1,protein2)), collapse = '_'))
string.full.use <- string.full[string.type$pair]
string.full.use.unique <- string.full.use[,3:17] %>% unique()

setDT(string.full)
setkey(string.full,pair)


string.detailed.use <- string.detailed[string.type$pair]



string.detailed.use.unique <- string.detailed.use[,3:11] %>%
     unique()
length(unique(string.detailed.use.unique$pair))
string.type.unique <- select(string.type, c('pair','score')) %>%
     unique() %>%
     left_join(string.detailed.use.unique, by = 'pair') %>%
     arrange(pair)
length(unique(string.type.unique$pair))

View(filter(string.detailed.use, cooccurence > 0))

colns <- colnames(string.detailed)[3:9]
prior <- 0.041

pr <- '9606.ENSP00000251810_9606.ENSP00000300738'
pr <- '9606.ENSP00000001146_9606.ENSP00000360958'
pr <- '9606.ENSP00000263025_9606.ENSP00000282908'
pr <- '9606.ENSP00000000233_9606.ENSP00000256682'
pr <- '9606.ENSP00000001008_9606.ENSP00000221114'
pr <- '9606.ENSP00000000233_9606.ENSP00000222547'
pr <- '9606.ENSP00000001146_9606.ENSP00000333212'

pr <- '9606.ENSP00000000233_9606.ENSP00000248901'
View(filter(string.type, pair == pr))
View(filter(string.detailed.use, pair == pr))
View(filter(string.full.use, pair == pr))
# View(filter(string.full, pair == pr))

calc.score <- function(df){
     df1 <- filter(string.detailed.use, pair == pr)[,3:9]
     #df1$cooccurence <- 0
     df <- as.matrix(df1)
     #df2 <- matrix(rep(c(0,0,0,621,585,0,0),2), ncol = 7, byrow = T)
     mat <- (df1/1000-prior)/(1-prior)
     mat[mat < 0] <- 0
     (1-exp(rowSums(log(1-mat))))+prior*exp(rowSums(log(1-mat)))
}


