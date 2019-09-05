
############################################################################

################## DICE

fls <- list.files('Raw/DICE/')
fls <- fls[-grep('merged', fls)]
fls

cells <- c('B cell, naive', 'T cell, CD4, naive', 'T cell, CD4, naive [activated]',
           'T cell, CD8, naive', 'T cell, CD8, naive [activated]', 'Monocyte, non-classical',
           'Monocyte, classical', 'NK cell, CD56dim CD16+', 'T cell, CD4, TFH', 
           'T cell, CD4, TH1', 'T cell, CD4, TH17', 'T cell, CD4, TH2',
           'T cell, CD4, TH1/17', 'T cell, CD4, memory TREG', 'T cell, CD4, naive TREG')

tpmList <- list()
CellType <- c()

for (i in 1:15) {
  flName <- fls[i]
  sam <- gsub('_TPM.csv', '', flName)
  print (sam)
  
  samples <- c(samples, sam)
  fl <- file.path('Raw/DICE/', flName)
  tpm <- read.csv(fl, header = T, stringsAsFactors = F)
  
  if (i==1) {
    refGenes <- sapply(tpm$Additional_annotations, function(x) strsplit(x, ';', fixed=T)[[1]][1])
  }
  
  genes <- sapply(tpm$Additional_annotations, function(x) strsplit(x, ';', fixed=T)[[1]][1])
  print(sum(genes != refGenes))
  
  tpm <- tpm[,-c(1:3)]
  
  filter <- which(duplicated(genes))
  tpm <- tpm[-filter,]
  genes <- genes[-filter]
  
  rownames(tpm) <- genes
  colnames(tpm) <- paste(sam, 1:ncol(tpm), sep='_')
  tpmList[[cells[i]]] <- tpm
  
  CellType <- c(CellType, rep(cells[i], ncol(tpm)))
  
}

tpmMatrix <- do.call(cbind, tpmList)

dim(tpmMatrix)
tpmMatrix[1:5,1:5]

CellType <- sapply(colnames(tpmMatrix), function(x) strsplit(x, '.', fixed=T)[[1]][1])
CellType

colnames(tpmMatrix) <- sapply(colnames(tpmMatrix), function(x) strsplit(x, '.', fixed=T)[[1]][2])
tpmMatrix <- round(tpmMatrix, 2)
saveRDS(tpmMatrix, file='rData/DICE_TPM.rds')

exprData <- tpmMatrix
phenoData <- data.frame(Database='DICE', 
                        CellType=CellType, 
                        row.names=colnames(exprData))
phenoData

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file='rData/eSet_DICE_TPM.rds')




############################################################################

################## BLUEPRINT

metadata <- read.table('Raw/Blueprint/BLUEPRINT_RNAseq_metadata.txt', header=T, sep='\t',
                       stringsAsFactors = F, fill = TRUE)
                              
unique(metadata$Sub.group)
unique(metadata$Cell.type)

fls <- file.path('Raw/Blueprint/RNAseq', metadata$fileName, 
                 fsep = .Platform$file.sep)

#expr <- read.table(filenames[1], header=T, stringsAsFactors = F)
#expr[1:5,1:5]
#expr$TPM

#tpmMatrix <- do.call("cbind", lapply(fls, function(fl) 
#  read.table(fl, header=T, stringsAsFactors = F)$TPM))
#tpmMatrix[1:5,1:5]
#rownames(tpmMatrix) <- read.table(fls[1], header=T, stringsAsFactors = F)$gene_id
#rownames(tpmMatrix) <- unlist(lapply(strsplit(rownames(tpmMatrix), 
#                                              '.', fixed=TRUE), function(gene) gene[1]))
#colnames(tpmMatrix) <- metadata$sampleName

tpmMatrix <- c()
for (i in 1:length(fls)) {
  flName <- fls[i]
  #print (flName)
  tpm <- read.table(flName, header = T, stringsAsFactors = F)
  
  if (i==1) {
    refGenes <- sapply(tpm$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
  }
  
  genes <- sapply(tpm$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
  print(sum(genes != refGenes))
  tpmMatrix <- cbind(tpmMatrix, tpm$TPM)
}

rownames(tpmMatrix) <- genes
colnames(tpmMatrix) <- metadata$sampleName[1:5]
saveRDS(tpmMatrix, 'rData/BLUEPRINT_TPM_ENSEMBL.rds')


#########
phenoData <- read.table('Raw/Blueprint/20160816.data.indexV2.txt', header=T, stringsAsFactors = F,
                        sep='\t', na.strings = '-')
phenoData[1:5,]

rownames(phenoData) <- phenoData$SAMPLE_NAME
phenoData <-phenoData[colnames(tpmMatrix),]
rownames(phenoData) == colnames(tpmMatrix)
                  
#phenoData$CELL_TYPE[which(is.na(phenoData$CELL_TYPE))] <-phenoData$BIOMATERIAL_TYPE[which(is.na(phenoData$CELL_TYPE))]
#phenoData$CELL_TYPE <- tolower(phenoData$CELL_TYPE)

View(phenoData)

exprData <- tpmMatrix
dim(exprData)

annoData <- readRDS(file='~/Projects/AnnotationData/annotation.gencode.v29.rds')
annoData

ovlp <- intersect(rownames(exprData), rownames(annoData))
ovlp

exprData <- exprData[ovlp,]
annoData <- annoData[ovlp,]

filter <- which(duplicated(annoData$gene_name))
filter

exprData <- exprData[-filter,]

rownames(exprData) <- annoData[-filter,]$gene_name
dim(exprData)
saveRDS(exprData, 'rData/BLUEPRINT_TPM_SYMBOL.rds')

eSet <- ExpressionSet(assayData = as.matrix(exprData),
                      phenoData = AnnotatedDataFrame(phenoData))

saveRDS(eSet, file='rData/BLUEPRINT_TPM_SYMBOL_eSet.rds')
