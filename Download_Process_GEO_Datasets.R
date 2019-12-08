###################### Download and process datasets
library(GEOquery)
library(oligo)

gse <- 'GSE46691'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE
#seriesMatrix <- getGEO(filename="~/Downloads/GSE60341_series_matrix.txt.gz")

phenoData <- pData(seriesMatrix[[1]])
phenoData

phenoData <- phenoData[,c(1,2,8,36:38)]
colnames(phenoData)

colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
colnames(phenoData) <- gsub(' ', '_', colnames(phenoData))

colnames(phenoData)[3] <- 'tissue'

seriesMatrix[[1]]@featureData@data

View(phenoData)

exprData <- exprs(seriesMatrix[[1]])
View(exprData)

# Ensembl
# CDF: O
install.packages("http://mbni.org/customcdf/24.0.0/ensg.download/pd.huex10st.hs.ensg_24.0.0.tar.gz",
                 repos = NULL, type = "source")

# Gencode
install.packages("http://mbni.org/customcdf/24.0.0/gencodeg.download/pd.huex10st.hs.gencodeg_24.0.0.tar.gz",
                 repos = NULL, type = "source")

library(pd.huex10st.hs.ensg)
library(pd.huex10st.hs.gencodeg)

filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE,filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg') #pd.huex.1.0.st.v2

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
