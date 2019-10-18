library(stringr)
library(rtracklayer)
library(biomaRt)

############# RefSeq (ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/)
gtf <- readGFF('data/ref_GRCm38_top_level.gff3.gz') # For featureCounts mm10
gtf <- gtf[with(gtf, type=='gene'),]

entrez <- unlist(lapply(gtf$Dbxref, function(x) paste(x, collapse = ';')))
entrez <- gsub('GeneID:', '', str_extract(entrez, 'GeneID:\\d+'))

annoRefSeq <- data.frame(Symbol=gtf$Name, EntrezID=as.character(entrez),
                         stringsAsFactors = F)

############# biomaRt

listMarts()

ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)

ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
attributes <- ensembl@attributes

annoBiomaRt <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', #'hgnc_symbol', 
                                  #'external_gene_name', 
                                  'mgi_symbol', 'description', 'gene_biotype'),
                     #filters = 'affy_hg_u133_plus_2', 
                     #values = affyids, 
                     mart = ensembl)


############ MGI (http://www.informatics.jax.org/downloads/reports/index.html)

#### MGI query
annoMGI <- read.table('data/MGIBatchReport_20191017_014740.txt',
                      header = T, sep = '\t', stringsAsFactors = F)
annoMGI

#### RPT from MGI (http://www.informatics.jax.org/downloads/reports/index.html)
annoRPT <- read.table('data/MGI_Gene_Model_Coord.rpt',
                      header = T, sep = '\t', stringsAsFactors = F)

colnames(annoRPT)[1:14] <- colnames(annoRPT)[2:15]
colnames(annoRPT) <- gsub('X\\d+..', '', colnames(annoRPT))

annoRPT <- annoRPT[,-15]

#### GFF3 from MGI (http://www.informatics.jax.org/downloads/reports/index.html)

gtf <- readGFF('data/E-MTAB-7462/MGI.gff3')
gtf <- gtf[with(gtf, type=='gene'),]

entrez <- unlist(lapply(gtf$Dbxref, function(x) paste(x, collapse = ';')))
entrez
entrez <- gsub('NCBI_Gene:', '', str_extract(entrez, 'NCBI_Gene:\\d+'))

annoGFF <- data.frame(Symbol=gtf$Name, EntrezID=as.character(entrez),
                      stringsAsFactors = F)
