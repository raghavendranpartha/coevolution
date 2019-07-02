require(VennDiagram)

hsap.genes <- readLines('../data/hsap.genes')
hsap.vert.genes <- readLines('../data/hsap.vert.genes')
hsap.dmel.genes <- readLines('../data/hsap.dmel.genes')
hsap.cele.genes <- readLines('../data/hsap.cele.genes')
hsap.scer.genes <- readLines('../data/hsap.scer.genes')

vd <- venn.diagram(list(vert = hsap.vert.genes,
                        fly = hsap.dmel.genes,
                        yeast = hsap.scer.genes,
                        worm = hsap.cele.genes), NULL)
grid.draw(vd)


hg19.uniprot.wdmelortholog <- readLines('../data/inparanoid.ortholog.maps/dmel_hsap/cleanedparsedOutput.D.melanogaster-H.sapiens.hsap.uniprot.id')
hg19.uniprot.wscerortholog <- readLines('../data/inparanoid.ortholog.maps/hsap_scer/cleanedparsedOutput.H.sapiens-S.cerevisiae.hsap.uniprot.id')
hg19.uniprot.wceleortholog <- readLines('../data/inparanoid.ortholog.maps/cele_hsap/cleanedparsedOutput.C.elegans-H.sapiens.hsap.uniprot.id')

vd2 <- venn.diagram(list(dmel = hg19.uniprot.wdmelortholog,
                        scer = hg19.uniprot.wscerortholog,
                        cele = hg19.uniprot.wceleortholog), NULL)
grid.draw(vd2)
