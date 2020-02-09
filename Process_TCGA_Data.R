
setwd('~/bigdata/TCGA/')

BiocManager::install('GDCRNATools')

library(GDCRNATools)
library(rtracklayer)
library(tibble)

#===============================================================================================
project <- 'TCGA-PRAD'
rnadir <- paste('data', project, 'RNAseq', sep='/')
#mirdir <- paste('data', project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)


####### Download mature miRNA data #######
#gdcRNADownload(project.id     = project, 
#               data.type      = 'miRNAs', 
#               write.manifest = FALSE,
#               method         = 'gdc-client',
#               directory      = mirdir)



####### Download clinical data #######
clinicaldir <- paste('data', project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = project, 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)



####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = project,
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)



####### Parse miRNAs metadata #######
#metaMatrix.MIR <- gdcParseMetadata(project.id = project,
#                                   data.type  = 'miRNAs', 
#                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
#metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
#metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

View(rnaCounts[1:5,])

saveRDS(rnaCounts, file=paste0('data/RNAseq_Counts_', gsub('-', '_', project), '.RDS'))


####### Merge miRNAs data #######
#mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
#                         path      = mirdir, # the folder in which the data stored
#                         organized = FALSE, # if the data are in separate folders
#                         data.type = 'miRNAs')

#saveRDS(rnaCounts, file=paste0('data/miRNAs_Counts_', gsub('-', '_', project), '.RDS'))



####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

View(clinicalDa)

saveRDS(clinicalDa, file=paste0('data/Clinical_', gsub('-', '_', project), '.RDS'))



#########################################################################################


## submitter id

gdcRNAMerge <- function(metadata, path, data.type, organized=FALSE) {
  
  #if (endsWith(path, '/')) {
  #  path = substr(path, 1, nchar(path)-1)
  #}
  
  if (organized==TRUE) {
    filenames <- file.path(path, metadata$file_name, 
                           fsep = .Platform$file.sep)
  } else {
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
  }
  
  if (data.type=='RNAseq') {
    message ('############### Merging RNAseq data ################\n',
             '### This step may take a few minutes ###\n')
    
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.table(gzfile(fl))$V2))
    rownames(rnaMatrix) <- read.table(gzfile(filenames[1]))$V1
    rownames(rnaMatrix) <- unlist(lapply(strsplit(rownames(rnaMatrix), 
                                                  '.', fixed=TRUE), function(gene) gene[1]))
    #colnames(rnaMatrix) <- metadata$sample
    colnames(rnaMatrix) <- metadata$submitter_id
    
    #rnaMatrix <- rnaMatrix[biotype$ensemblID,]
    
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
  } else if (data.type=='pre-miRNAs') {
    message ('############### Merging pre-miRNAs data ################\n',
             '### This step may take a few minutes ###\n')
    
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.delim(fl)$read_count))
    rownames(rnaMatrix) <- read.delim(filenames[1])$miRNA_ID
    
    colnames(rnaMatrix) <- metadata$sample
    
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
    
    
  } else if (data.type=='miRNAs') {
    message ('############### Merging miRNAs data ###############\n')
    
    mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
    #mirs <- sort(unique(names(unlist(mirMatrix))))
    mirs <- rownames(mirbase)
    mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                         function(expr) expr[mirs]))
    
    rownames(mirMatrix) <- mirbase$v21[match(mirs,rownames(mirbase))]
    colnames(mirMatrix) <- metadata$sample
    
    mirMatrix[is.na(mirMatrix)] <- 0
    
    nSamples = ncol(mirMatrix)
    nGenes = nrow(mirMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of miRNAs: ', nGenes, '\n', sep=''))
    
    return (mirMatrix)
  } else {
    return ('error !!!')
  }
}


