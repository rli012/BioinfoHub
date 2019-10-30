immune.signature1 <- readRDS(file='data/immune.signature1.RDS')
immune.signature1$Genes <- gsub('/IL8|@', '', immune.signature1$Genes)

immune.signature1$Genes <- gsub(',\\s*', '|', toupper(immune.signature1$Genes))
View(immune.signature1)

signatures <- sapply(immune.signature1$Genes, function(x) sort(unique(strsplit(x, '|', fixed=T)[[1]])))

names(signatures) <- immune.signature1$Subcluster
signatures

summary(gsvaData)



gsvaData <- gsva(expr = exprLogCPM,
                 gset.idx.list = signatures,
                 method="gsva")


fit <- lmFit(gsvaData, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
