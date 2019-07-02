library(RERconverge)
require(dplyr)
require(data.table)
hsap.genes <- readLines('../data/hsap.genes')
hsap.vert.genes <- readLines('../data/hsap.vert.genes')
hsap.dmel.genes <- readLines('../data/hsap.dmel.genes')
hsap.cele.genes <- readLines('../data/hsap.cele.genes')
hsap.scer.genes <- readLines('../data/hsap.scer.genes')

omaage <- 