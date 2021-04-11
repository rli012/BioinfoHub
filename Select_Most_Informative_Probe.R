
### select most informative probe, MAX IQR
library(plyr)

iqr <- apply(expr[,-ncol(expr)], 1, IQR)
iqr

expr$IQR <- iqr


probeIndex <- ddply(expr, .(ID_REF), summarise, probe=which.max(IQR))
probeIndex[1:5,]

probes <- c()
for (i in 1:nrow(probeIndex)) {
  probe <- rownames(expr)[which(expr$ID_REF==probeIndex$ID_REF[i])]
  print (probe)

  probes <- c(probes, probe[probeIndex$probe[i]])
}

filter <- which(is.na(probes))
filter

#probes <- probes[-filter]
#probes

exprData <- exprData[probes,]
exprData
