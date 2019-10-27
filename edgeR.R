### edgeR

dge <-  DGEList(counts = countMatrix)
rpkmData <- rpkm(dge, gene.length = geneLength)

#keep <- rowSums(rpkmData > 1) >= 2 #0.5*ncol(countMatrix)
#dge <- dge[keep,,keep.lib.sizes = FALSE]

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 8 #0.5*ncol(countMatrix)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = FALSE]


comparison <- 'sAg_loss_BASELINE - No_sAg_loss_BASELINE'

contrast.matrix <- makeContrasts(contrasts=comparison,
                                 levels=design)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, contrast=contrast.matrix)

dgeTable <- topTags(qlf, n = Inf, sort.by = 'p', adjust.method = 'BH')
dgeTable <- dgeTable[[1]]
dgeTable$Symbol <- rownames(dgeTable)

View(dgeTable)
