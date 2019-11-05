dge <-  DGEList(counts = countMatrix)
rpkmData <- rpkm(dge, gene.length = geneLength)

keep <- rowSums(edgeR::cpm(dge) > 1) >= 8 #0.5*ncol(countMatrix)
sum(keep)

rpkmData <- rpkmData[keep,]

write.table(rpkmData, file='report/RPKM_for_CYBERSORT.txt', quote=F, sep='\t')
write.table(rpkmData, file='report/RPKM_for_ssGSEA.txt', quote=F, sep='\t')
