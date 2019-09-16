
library(Biobase)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)


setwd('~/Projects/Infrastructure_20190109/BioinfoHub/')

################################ HELP FUNCTIONS

organizeEnrichFun <- function(go) {

  Terms <- paste(go$ID, go$Description, sep='~')
  Counts <- go$Count

  GeneRatio <- go$GeneRatio
  BgRatio <- go$BgRatio

  pValue <- go$pvalue
  FDR <- go$p.adjust

  listTotal <- vapply(go$GeneRatio, function(v)
    convertRatioFun(v, type='bg'), numeric(1))
  popHits <- vapply(go$BgRatio, function(v)
    convertRatioFun(v, type='hit'), numeric(1))
  popTotal <- vapply(go$BgRatio, function(v)
    convertRatioFun(v, type='bg'), numeric(1))

  foldEnrichment <- as.vector(Counts/listTotal*popTotal/popHits)

  geneID <- go$geneID
  #geneSymbol <- unlist(lapply(strsplit(geneID, '/', fixed=TRUE),
  #                            function(v) paste(ensembl2symbolFun(v), collapse = '/')))

  goOutput <- data.frame(Terms, Counts, GeneRatio, BgRatio, pValue, FDR,
                         foldEnrichment, geneID)#, geneSymbol)

  return (goOutput)
}


###
convertRatioFun <- function(v, type='bg') {
  ratio <- strsplit(v, '/', fixed=TRUE)

  if (type=='bg') {
    num <- as.numeric(as.character(ratio[[1]][2]))
  } else if (type=='hit') {
    num <- as.numeric(as.character(ratio[[1]][1]))
  }

  return (num)
}


#####################################################################

########## Data preprocessing and Differential gene expression

### Import the ExpressionSet
eSet <- readRDS('data/TCGA_RNAseq_toy_eSet.rds')

#eSet <- readRDS('data/GSE79209_Raw_eSet.rds')

### Counts data
countsMatrix <- exprs(eSet)
countsMatrix[1:5,1:5]

### Phenotypic data
phenoData <- pData(eSet)
View(phenoData)

### Annotation data
annoData <- eSet@featureData@data

### Create DGEList object
dge <-  DGEList(counts = countsMatrix)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(cpm(dge) > 1) >= 0.5*ncol(countsMatrix)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
v <- voom(dge, design=NULL, plot = FALSE)

exprAfterVoom <- v$E ### for visualization
exprLogCPM <- cpm(dge,log = TRUE) ### for visualization
exprLogCPM

### Prepare comparison matrix
group <- factor(phenoData$sample_type)
#group <- factor(phenoData$HistologyGrade)

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(contrasts='PrimaryTumor - SolidTissueNormal',
                                 levels=design)
contrast.matrix

### Differential gene expression analysis (limma)

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


### Report DEGs
dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
#dgeTable <- topTable(fit2, coef='Moderate', n=Inf, adjust.method='BH', sort.by='p')
View(dgeTable)


### Map Ensembl ID to gene symbol
idx <- match(rownames(dgeTable), rownames(annoData))
dgeTable$Symbol <- annoData$Symbol[idx]
dgeTable$Biotype <- annoData$Biotype[idx]


###############################################################################

########## Volcano Plot

adjPvalThreshold <- 0.05
logFcThreshold <- log2(2)

dataForVolcanoPlot <- dgeTable

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, abs(logFC) < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'

#dataForVolcanoPlot$Symbol <- toupper(dataForVolcanoPlot$Symbol)


g <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color=Significance), alpha=1, size=1) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), color='darkgreen',linetype='dashed')+
  scale_color_manual(values = c("green3", "black", "red")) +
  xlab(expression('log'[2]*'(Fold Change)')) + ylab(expression('-log'[10]*'(P Value)')) +
  #xlim(-11,11) + #ylim(0,2.5) +
  geom_text_repel(data = subset(dataForVolcanoPlot, -log10(adj.P.Val)>10 & logFC>0),
                  segment.alpha = 0.4, color='red',
                  aes(label = Symbol), size = 3, segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot, -log10(adj.P.Val)>10 & logFC<0),
                  segment.alpha = 0.4, color='green3',
                  aes(label = Symbol), size = 3, segment.color = 'black') +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=16))

print(g)

############## BarPlot

up <- sum(dataForVolcanoPlot$Significance =='UP')
down <- sum(dataForVolcanoPlot$Significance =='DOWN')

dataForBarPlot <- data.frame(geneNums = c(up, down),
                             Regulation = factor(c('Up-regulated','Down-regulated'),
                                                 levels=c('Up-regulated','Down-regulated')))


ggplot(data=dataForBarPlot, aes(x=Regulation, y=geneNums,
                   fill=Regulation)) + geom_bar(stat = 'identity') +
  scale_x_discrete(limits=dataForBarPlot$geneClass) +
  scale_fill_discrete(name = "") +
  ylab('No. of Differentially Expressed Genes') +xlab('') +
  #theme(axis.text.x = element_text(angle = angle)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='white'),
                   panel.background = element_blank(),
                   legend.position='none',
                   axis.title.y = element_text(size=16),
                   axis.text.x = element_text(angle = 0,size=14),
                   axis.text.y=element_text(size=14))


########### Boxplot

gene <- 'INSIG1'
ensembl <- rownames(annoData)[which(annoData$Symbol==gene)]
ensembl

expr <- exprLogCPM[ensembl,]
expr

group <- factor(phenoData$sample_type, levels=c('SolidTissueNormal', 'PrimaryTumor'))
group

dataForBoxPlot <- data.frame(expr, group)
dataForBoxPlot

ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_boxplot(aes(fill=group),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_smooth(method='lm') +
  #facet_wrap(~dataset) +
  geom_jitter(size=2, width=0.05, color='blue') +
  xlab('') + ylab(expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16))





########### Heatmap

ann_colors <- list(
  Sample_Type=c(PrimaryTumor='lightcoral',
              SolidTissueNormal='darkolivegreen')
)

anno = HeatmapAnnotation(Sample_Type=phenoData$sample_type,
                          col=ann_colors,
                          simple_anno_size_adjust = TRUE,
                          #annotation_height = c(1,1),
                          height = unit(8, "mm"),
                          #summary = anno_summary(height = unit(4, "cm")),
                          show_legend = c("bar" = TRUE),
                          show_annotation_name = FALSE)


genes <- rownames(dataForVolcanoPlot)[dataForVolcanoPlot$Significance != 'NS']
genes

scaledExprData <- t(scale(t(exprLogCPM[genes,])))
scaledExprData[1:5,1:5]


ht <- Heatmap(scaledExprData, name = 'Expression',
              top_annotation = anno,
              #column_split = factor(phenoData$sample_type,
              #                      levels=str_sort(unique(phenoData$sample_type), numeric = T)),
              column_title = NULL,
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              show_row_dend = FALSE,
              column_names_rot = 90,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 12))

draw(ht)



#####################################################################

########### GO, KEGG Pathway Analysis
### clusterProfiler (R/Bioconductor) or Metascape (http://metascape.org/gp/index.html#/main/step1)

genes <- rownames(dataForVolcanoPlot)[dataForVolcanoPlot$Significance != 'NS']
genes

dataForEnrichment <- annoData[genes,]
dataForEnrichment

# SLOW !!!
go <- list()
for (category in c('BP','CC','MF')) {
  go[[category]] <- enrichGO(gene = rownames(dataForEnrichment),
                             universe = rownames(annoData),
                             OrgDb = org.Hs.eg.db,
                             ont = category,
                             keyType = 'ENSEMBL',
                             pAdjustMethod = "fdr",
                             pvalueCutoff = 0.01,
                             readable = FALSE)

  #go[[category]] <- simplify(go[[category]], cutoff=0.7, by="p.adjust", select_fun=min)
  go[[category]] <- organizeEnrichFun(data.frame(go[[category]]@result))
  go[[category]]$Category <- category

}

goForPlot <- do.call(rbind, go)
goForPlot$Terms

goForPlot <- goForPlot[goForPlot$pValue<0.05,]
dim(goForPlot)

ggplot(data=goForPlot, mapping=aes(x=Terms, y=-log(pValue,10), fill=Category)) +
  geom_bar(stat='identity') + scale_x_discrete(limits=goForPlot$Terms) +
  ylim(0, max(-log(goForPlot$pValue,10))) + theme(legend.title=element_blank())+
  ylab(expression('-log'[10]*'(P Value)'))+ xlab('') + #coord_flip() +
  #scale_fill_hue(breaks=go$Category,labels=go$Category) + #name = "Domain",
  scale_fill_manual(name='', values=c('darkblue','darkgreen','darkred'),
                    labels=c('Biological Process','Cellular Component','Molecular Function'))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='white'),
                   panel.background = element_blank())+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        axis.text.x=element_text(angle=45, hjust=1)) +
  theme(legend.text = element_text(size = 10),
        legend.position = 'top')+
  theme(plot.margin=margin(t = 0.25, r = 0.25, b = 0.25, l = 6, unit = "cm"))


###

kegg <- enrichKEGG(gene = as.character(unique(dataForEnrichment$Entrez)),
                   organism = 'hsa',
                   universe = as.character(unique(annoData$Entrez)),
                     #dataForEnrichment$Entrez[!is.na(dataForEnrichment$Entrez)])),
                   minGSSize = 5,
                   maxGSSize = 500,
                   pAdjustMethod = 'fdr',
                   pvalueCutoff = 0.01)

kegg <- data.frame(kegg@result)
kegg

kegg$geneID <- unlist(lapply(kegg$geneID, function(v)
  paste(rownames(dataForEnrichment)[match(strsplit(v, '/', fixed=TRUE)[[1]],
                                       dataForEnrichment$Entrez)], collapse = '/')))

kegg$Category='KEGG'

keggForPlot <- kegg[kegg$pvalue<0.05,]
keggForPlot$Terms <- paste(keggForPlot$ID, keggForPlot$Description, sep='~')



ggplot(data=keggForPlot, mapping=aes(x=Terms, y=-log(pvalue,10))) +
  geom_bar(stat='identity',fill='black') + scale_x_discrete(limits=rev(keggForPlot$Terms)) +
  ylim(0, max(-log(kegg$pvalue,10))) + theme(legend.title=element_blank())+ylab('-log10(P Value)')+
  xlab('') + coord_flip() +
  scale_fill_manual(values=c('black'))+#,breaks=kegg$Reg) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  #theme(legend.text = element_text(size = 12)) +
  theme(legend.position = 'none') +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='white'),
                   panel.background = element_blank())


