library(biomaRt)

listMarts()

ensembl=useMart("ensembl")
ensembl
datasets <- listDatasets(ensembl)
head(datasets)

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- ensembl@attributes
View(attributes)

#affyids=c("202763_at","209310_s_at","207500_at")
annotation <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', #'hgnc_symbol', 
                                 'external_gene_name', 'description', 'gene_biotype'),
                    #filters = 'affy_hg_u133_plus_2', 
                    #values = affyids, 
                    mart = ensembl)

annotation <- getBM(attributes=c('ensembl_gene_id', 'uniprot_gn_id', 'uniprotswissprot'),
                    #filters = 'affy_hg_u133_plus_2', 
                    #values = affyids, 
                    mart = ensembl)

dim(annotation)
View(annotation)

searchDatasets(mart = ensembl, pattern = "hsapiens")
searchAttributes(mart = ensembl, pattern = "hgnc")

### Ensembl Archive
listEnsemblArchives()
listMarts(host = 'oct2016.archive.ensembl.org')

ensembl86 <- useMart(host='oct2016.archive.ensembl.org', 
                     biomart='ENSEMBL_MART_ENSEMBL', 
                     dataset='hsapiens_gene_ensembl')
