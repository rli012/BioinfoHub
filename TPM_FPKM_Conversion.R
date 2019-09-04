
### EXAMPLE

geneLengths <- readRDS('~/Projects/Infrastructure/data/gencode.v22.annotation.gene_lengths.RDS')

tcga <- read.table(gzfile("~/Projects/Infrastructure/data/7d9d7119-4a9e-46f6-b758-153bfcc869a5.htseq.counts.gz"),
                   header=F, stringsAsFactors = F, row.names = 1)

geneOvlp <- Reduce(intersect, list(rownames(geneLengths), rownames(tcga)))
geneOvlp


tcga$V1 <- tcga$V2*2
tcga <- tcga[geneOvlp,]
head(tcga)

sum(rownames(tcga) != rownames(geneLengths))

tcgaNorm <- tcga/geneLengths$Length*1e3
tcgaTPM <- sweep(tcgaNorm,2,colSums(tcgaNorm),`/`)*1e6
head(tcgaTPM)


tcgaNorm <- sweep(tcga,2,colSums(tcga),`/`)*1e6
tcgaFPKM <- tcgaNorm/geneLengths$Length*1e3
head(tcgaFPKM)


testFPKM <- read.table(gzfile("~/Projects/Infrastructure/data/7d9d7119-4a9e-46f6-b758-153bfcc869a5.FPKM.txt.gz"),
                       header=F, stringsAsFactors = F)
rownames(testFPKM) <- testFPKM$V1
head(testFPKM[geneOvlp,])



### only use reads mapped to protein coding genes
pc <- readRDS('~/Projects/Infrastructure/data/gencode.v22.protein_coding.RDS')
pc


tcgaNorm <- sweep(tcga,2,colSums(tcga[pc,], na.rm=T),`/`)*1e6
tcgaFPKM <- tcgaNorm/geneLengths$Length*1e3
head(tcgaFPKM)

head(testFPKM[geneOvlp,])

sum(round(tcgaFPKM$V2,6) == round(testFPKM[geneOvlp,]$V2,6))
sum(round(tcgaFPKM$V2,6) != round(testFPKM[geneOvlp,]$V2,6))


tcgaNorm <- tcga/geneLengths$Length*1e3
tcgaTPM <- sweep(tcgaNorm,2,colSums(tcgaNorm[pc,]),`/`)*1e6
head(tcgaTPM)



####################################################################

###### GENERATE GENE LENGTHS

library(GenomicFeatures)

gtf_file <- "~/Projects/Infrastructure/data/gencode.v22.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

saveDb(txdb, file="~/Projects/Infrastructure/data/gencode.v22.annotation.sqlite")

txdb <- loadDb("~/Projects/Infrastructure/data/gencode.v22.annotation.sqlite")
txdb


tcga <- read.table(gzfile("~/Projects/Infrastructure/data/7d9d7119-4a9e-46f6-b758-153bfcc869a5.htseq.counts.gz"),
                   header=F, stringsAsFactors = F, row.names = 1)
#rownames(tcga) <- sapply(rownames(tcga), function(g) strsplit(g, '.', fixed=T)[[1]][1])
tcga$V1 <- tcga$V2*2
head(tcga)


vennData <- list(GENCODE22=names(genes(txdb)), TCGA=rownames(tcga))
vennDiagramColors <- c('red', 'blue')#, 'green', 'black', 'orange')#, 'purple')
venn.diagram(vennData, filename = '~/Projects/Infrastructure/data/Venn.TCGA.GENCODE22.png', 
             imagetype = 'png', col = vennDiagramColors, lwd = 2, fill = vennDiagramColors, alpha = rep(0.25, length(vennDiagramColors)),
             fontfamily = 'sans', margin = 0.2, cat.dist = 0.03, cat.col = vennDiagramColors, cat.fontface = 'bold',
             cat.fontfamily = 'sans', cat.cex = 0.8, cex = 0.75)


geneOvlp <- Reduce(intersect, list(names(genes(txdb)), rownames(tcga)))
length(geneOvlp)


geneDiff <- setdiff(rownames(tcga), geneOvlp)
geneDiff

#gr <- genes(txdb)[geneOvlp]
#gr

exons.list.per.gene <- exonsBy(txdb,by="gene")
exons.list.per.gene

# long non coding example
exons.list.per.gene[['ENSG00000243485.3']]


### SLOW !!!
exonic.gene.sizes <- mclapply(exons.list.per.gene,function(x){sum(width(reduce(x)))}, mc.cores=4)
exonic.gene.sizes

geneLengths <- data.frame(GeneID=names(exonic.gene.sizes), Length=unlist((exonic.gene.sizes)))
head(geneLengths)

geneLengths <- geneLengths[geneOvlp,]


### FAST, BUT NEED A .BAM FILE TO RUN
library(Rsubread)
fc <- featureCounts(files='wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam',
                    annot.ext='gencode.v22.annotation.gtf',
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType='exon',
                    GTF.attrType='gene_id')


geneLengths <- data.frame(GeneID=fc$annotation$GeneID, Length=fc$annotation$Length, row.names=fc$annotation$GeneID)
geneLengths

geneLengths <- geneLengths[geneOvlp,]
geneLengths[1:10,]


saveRDS(geneLengths, file = '~/Projects/Infrastructure/data/gencode.v22.annotation.gene_lengths.RDS')


tcga <- tcga[geneOvlp,]
tcga


tcgaNorm <- tcga/geneLengths$Length*1e3
tcgaTPM <- sweep(tcgaNorm,2,colSums(tcgaNorm),`/`)*1e6
head(tcgaTPM)



tcgaNorm <- sweep(tcga,2,colSums(tcga),`/`)*1e6
tcgaFPKM <- tcgaNorm/geneLengths$Length*1e3
head(tcgaFPKM)



testFPKM <- read.table(gzfile("~/Projects/Infrastructure/data/7d9d7119-4a9e-46f6-b758-153bfcc869a5.FPKM.txt.gz"),
                       header=F, stringsAsFactors = F)
rownames(testFPKM) <- testFPKM$V1
head(testFPKM[geneOvlp,])



###

require(rtracklayer)
#gtf <- readGFF(gtf_file, version=2L)
#gtf
#saveRDS(gtf, 'gencode.v22.annotation.RDS')

gtf <- readRDS('~/Projects/Infrastructure/data/gencode.v22.annotation.RDS')
head(gtf)
dim(gtf)

gtf <- gtf[with(gtf, type=='gene' & gene_type=='protein_coding'),]
head(gtf)

proteinCodingGenes <- gtf$gene_id

saveRDS(proteinCodingGenes, file='~/Projects/Infrastructure/data/gencode.v22.protein_coding.RDS')






########################################################################
###### ALTERNATIVE

#biocLite('scater')
library(scater)

#?calculateFPKM()
#?calculateTPM()

#biocLite('SingleCellExperiment')
library(SingleCellExperiment)


sce_tcga <- SingleCellExperiment(
  assays = list(counts = tcga), 
  colData = colnames(tcga)
)


tcga_tpm <- calculateTPM(object = sce_tcga, effective_length = geneLengths$Length)
tcga_tpm


