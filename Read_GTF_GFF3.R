library(rtracklayer)

gtf <- readGFF('gencode.v29.annotation.gtf.gz', version=2L)
head(gtf)
#saveRDS(gtf, 'annotation.gencode.v29.all.rds')

#gtf <- gtf[with(gtf, type=='gene' & gene_type=='protein_coding'),]
#head(gtf)
#proteinCodingGenes <- gtf$gene_id

gtf <- gtf[with(gtf, type=='gene'),]
head(gtf)
dim(gtf)

genes <- sapply(gtf$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
filter <- which(duplicated(genes))

gtf <- gtf[-filter,]
rownames(gtf) <- genes[-filter]

saveRDS(gtf, file='annotation.gencode.v29.rds')
