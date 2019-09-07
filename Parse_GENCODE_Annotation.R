
library(rtracklayer)
library(tibble)

getGENCODEAnnotation <- function(species='human', release='31', type='gene') {

  # species: human, mouse
  # release: human 31, mouse M20
  # type: gene, transcript
  
  baseurl <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
  
  gtf.file <- paste0(baseurl, 'Gencode_', species, '/release_', release, '/gencode.v', release, '.annotation.gtf.gz')
  gtf <- readGFF(gtf.file, version=2L)
  
  if (type!='all') {
    gtf <- gtf[gtf$type==type,]
    ensembl <- sapply(gtf$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
    gtf <- add_column(gtf, ensembl, .before = 'gene_id')
  }

  return(gtf)
}


gtf <- getGENCODEAnnotation(species='human', release='31', type='gene')
