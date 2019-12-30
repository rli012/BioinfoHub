
setwd('~/bigdata/PCa/data/fromSRA/GSE54460/')

library(Rsubread)

fls <- list.files('BAM/') # aligntest
fls <- fls[endsWith(fls, 'bam')]
fls <- paste0('BAM/', fls)
fls

### Ensembl
countMatrix <- featureCounts(files = fls, annot.ext = '~/bigdata/PCa/data/Reference/gencode.v32.annotation.gtf',
                             isGTFAnnotationFile = TRUE, GTF.attrType = 'gene_id',primaryOnly = TRUE, # minMQS = 255,
                             isPairedEnd = TRUE)

### RefSeq (In-built), include primary alignment of multiple mapping
countMatrix <- featureCounts(files = fls, annot.inbuilt = 'hg38', countMultiMappingReads = TRUE, 
                             primaryOnly = TRUE, isPairedEnd = TRUE)
