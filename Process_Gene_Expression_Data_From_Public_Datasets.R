
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

