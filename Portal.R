
rownames(phenoData) == colnames(countMatrix)

dim(countMatrix)
dim(phenoData)

ovlp <- intersect(rownames(phenoData),colnames(countMatrix))
ovlp

countMatrix <- countMatrix[,ovlp]
phenoData <- phenoData[ovlp,]

phenoData$PlotGroup <- paste0(phenoData$Disease, ' - ', phenoData$TissueSource)


### Create DGEList object
dge <-  DGEList(counts = countMatrix)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 5 #0.5*ncol(countMatrix)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = FALSE]

### Voom normalization
v <- voom(dge, design=NULL, plot = FALSE)

exprAfterVoom <- v$E ### for visualization
exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
exprLogCPM

exprData <- exprLogCPM

eSet <- new("ExpressionSet", exprs = exprData)
pData(eSet) <- phenoData

saveRDS(eSet, file=paste0('data/Portal/eSet.rds'))


