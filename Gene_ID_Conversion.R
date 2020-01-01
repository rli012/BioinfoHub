
library(rtracklayer)
library(tibble)
library(biomaRt)
library(stringr)
library(devtools)
#install_github("js229/Vennerable")
library(Vennerable)


#############################################################################################
# ==================================== Data Collection ==================================== #
#############################################################################################

####### GENCODE V32
getGENCODEAnnotation <- function(species='human', release='32', type='gene', gtf.file=NULL) {
  
  if (is.null(gtf.file)) {
    baseurl <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
    gtf.file <- paste0(baseurl, 'Gencode_', species, '/release_', release, '/gencode.v', release, '.annotation.gtf.gz')
  }
  
  gtf <- readGFF(gtf.file, version=2L)
  
  if (type!='all') {
    gtf <- gtf[gtf$type==type,]
    ensembl <- sapply(gtf$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
    gtf <- add_column(gtf, ensembl, .before = 'gene_id')
  }
  
  return(gtf)
}

# The gene annotation in "pseudoautosomal regions" (PAR) of chromosome Y is identical between chromosomes X and Y
# Gene names are usually HGNC or MGI-approved gene symbols mapped to the GENCODE genes by the Ensembl xref pipeline. 
# Sometimes, when there is no official gene symbol, the Havana clone-based name is used.
gencode <- getGENCODEAnnotation(species='human', release='32', type='gene')
rownames(gencode) <- ifelse(duplicated(gencode$ensembl), gencode$gene_id, gencode$ensembl) # PAR_Y, or filter out
gencode <- add_column(.data = gencode, .after = 'end', length=gencode$end-gencode$start+1)
saveRDS(gencode, file='data/Annotation/GENCODE_Annotation_Human_V32_20191230.RDS')

idx <- which(duplicated(gencode$gene_name))
View(gencode[idx,])
sort(unique(gencode[idx,]$gene_name)) # 113 genes



####### ENSEMBL V98
getENSEMBLAnnotation <- function(attributes=NULL) {
  
  if (is.null(attributes)) {
    attributes <- c('ensembl_gene_id','chromosome_name','start_position','end_position','strand',
                    'entrezgene_id','hgnc_id','hgnc_symbol','external_gene_name','description','gene_biotype')
  }
  
  #listMarts()
  
  ensembl=useMart("ensembl")
  datasets <- listDatasets(ensembl)

  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  #attributes <- ensembl@attributes
  #View(attributes)
  
  #affyids=c("202763_at","209310_s_at","207500_at")
  annotation <- getBM(attributes=attributes,
                      #filters = 'affy_hg_u133_plus_2', 
                      #values = affyids, 
                      mart = ensembl)
  
  return(annotation)
}

ensembl <- getENSEMBLAnnotation()
saveRDS(ensembl, file='data/Annotation/ENSEMBL_Annotation_Human_V98_20191230.RDS')



###
ensembl <- readRDS(file='data/Annotation/ENSEMBL_Annotation_Human_V98_20191230.RDS')

ensembl[ensembl==''] <- NA

keep <- which(ensembl$chromosome_name %in% c(1:22, 'MT', 'X', 'Y'))
ensembl <- ensembl[keep,]

#filter <- which(duplicated(ensembl$ensembl_gene_id))
#ensembl <- ensembl[-filter,]
#dim(ensembl)


####### HGNC

# wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -O hgnc_complete_set_20191230.txt
# wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/withdrawn.txt -O hgnc_withdrawn_20191230.txt

hgnc <- read.table(file='data/Annotation/hgnc_complete_set_20191230.txt', sep='\t', quote='', 
                   comment.char='', header = T, stringsAsFactors = F)

hgnc.withdrawn <- read.table(file='data/Annotation/hgnc_withdrawn_20191230.txt', sep='\t', quote='', 
                             comment.char='', header = T, stringsAsFactors = F)

hgnc[hgnc==''] <- NA
hgnc.withdrawn[hgnc.withdrawn==''] <- NA

hgnc.withdrawn$CURRENT_ID <- unlist(lapply(hgnc.withdrawn$MERGED_INTO_REPORT.S...i.e.HGNC_ID.SYMBOL.STATUS.,
                                           function(x) strsplit(x,'|', fixed=T)[[1]][1]))
hgnc.withdrawn$CURRENT_SYMBOL <- unlist(lapply(hgnc.withdrawn$MERGED_INTO_REPORT.S...i.e.HGNC_ID.SYMBOL.STATUS.,
                                           function(x) strsplit(x,'|', fixed=T)[[1]][2]))


####### NCBI
ncbi.gene.info <- read.table('data/Annotation/NCBI_Homo_sapiens.gene_info.20191230.gz', header = T, stringsAsFactors = F, 
                             sep='\t', quote='', comment.char = '')

ncbi.gene.info[ncbi.gene.info=='-'] <- NA

ncbi.gene.info$ensembl_gene_id <- str_extract(ncbi.gene.info$dbXrefs, 'ENSG\\d+')
ncbi.gene.info$hgnc_id <- str_extract(ncbi.gene.info$dbXrefs, 'HGNC:\\d+')

ncbi.gene.history <- read.table(file='data/Annotation/NCBI_Homo_sapiens.gene_history.20191230.txt', sep='\t', header=T,
                                stringsAsFactors = F, quote = '', comment.char = '')

ncbi.gene.history[ncbi.gene.history=='-'] <- NA
ncbi.gene.history$GeneID <- as.numeric(ncbi.gene.history$GeneID)
ncbi.gene.history$Discontinued_GeneID <- as.numeric(ncbi.gene.history$Discontinued_GeneID)





#############################################################################################
# ====================================== Update IDs ======================================= #
#############################################################################################

####### ENSEMBL

idx <- which(!ensembl$hgnc_id %in% c(hgnc$hgnc_id,NA))
idx
tmp <- hgnc.withdrawn[match(ensembl[idx,]$hgnc_id, hgnc.withdrawn$HGNC_ID),]
tmp

ensembl[idx,]$hgnc_id <- tmp$CURRENT_ID
ensembl[idx,]$hgnc_symbol <- tmp$CURRENT_SYMBOL
ensembl[idx,]$description <- NA

#which(ensembl$entrezgene_id %in% ncbi.gene.history$Discontinued_GeneID)
#which(!ensembl$entrezgene_id[which(!is.na(ensembl$entrezgene_id))] %in% ncbi.gene.info$GeneID)

idx <- which(ensembl$entrezgene_id %in% ncbi.gene.history$Discontinued_GeneID)
ensembl[idx,]$entrezgene_id

newids <- ncbi.gene.history$GeneID[match(ensembl[idx,]$entrezgene_id, ncbi.gene.history$Discontinued_GeneID)]
newids

ensembl[idx,]$entrezgene_id <- newids


####### HGNC

idx <- which(!hgnc$ensembl_gene_id %in% c(ensembl$ensembl_gene_id,NA))
idx

View(hgnc[idx,])

hgnc$ensembl_gene_id[idx] <- NA

which(hgnc$entrez_id[!is.na(hgnc$entrez_id)] %in% ncbi.gene.history$Discontinued_GeneID)
idx <- which(!hgnc$entrez_id %in% c(ncbi.gene.info$GeneID,NA))
idx

newids <- ncbi.gene.history$GeneID[match(hgnc[idx,]$entrez_id, ncbi.gene.history$Discontinued_GeneID)]
newids

hgnc[idx,]$entrez_id <- newids


####### NCBI

idx <- which(!ncbi.gene.info$ensembl_gene_id %in% c(ensembl$ensembl_gene_id,NA))
idx

dim(ncbi.gene.info)
sum(ncbi.gene.info$ensembl_gene_id %in% c(ensembl$ensembl_gene_id,NA))

ncbi.gene.info$ensembl_gene_id[idx] <- NA

#ensembl[which(ensembl$hgnc_id %in% hgnc.withdrawn$HGNC_ID),]

#ovlp <- intersect(hgnc$entrez_id, ncbi.gene.info$GeneID)
#ovlp

#gene1 <- hgnc$symbol[match(ovlp, hgnc$entrez_id)]
#gene2 <- ncbi.gene.info$Symbol[match(ovlp, ncbi.gene.info$GeneID)]

#data.frame(gene1, gene2)[which(gene1!=gene2),]

idx <- which(!ncbi.gene.info$hgnc_id %in% c(hgnc$hgnc_id,NA))
idx
ids <- hgnc.withdrawn[match(ncbi.gene.info[idx,]$hgnc_id, hgnc.withdrawn$HGNC_ID),]
ids

ncbi.gene.info[idx,]$hgnc_id <- ids$CURRENT_ID
ncbi.gene.info[idx,]$Symbol <- ids$CURRENT_SYMBOL
ncbi.gene.info[idx,]$description <- NA



##################

set1 <- ensembl$ensembl_gene_id[which(!is.na(ensembl$hgnc_id))]
set2 <- hgnc$ensembl_gene_id[which(!is.na(hgnc$ensembl_gene_id))]

set1 <- ensembl$ensembl_gene_id[which(!is.na(ensembl$entrezgene_id))]
set2 <- ncbi.gene.info$ensembl_gene_id[which(!is.na(ncbi.gene.info$ensembl_gene_id))]

set1 <- hgnc$entrez_id[which(!is.na(hgnc$entrez_id))]
set2 <- ncbi.gene.info$GeneID[which(!is.na(ncbi.gene.info$hgnc_id))]


set.list <- list(set1, set2)
vennData <- Venn(set.list)
vennData



#############################################################################################
# ======================================= Map IDs ========================================= #
#############################################################################################

####### ENSEMBL

idx<- which(is.na(ensembl$hgnc_id))
ensembl$hgnc_id[idx] <- hgnc$hgnc_id[match(ensembl$ensembl_gene_id[idx], hgnc$ensembl_gene_id)]

idx<- which(is.na(ensembl$entrezgene_id))
ensembl$entrezgene_id[idx] <- ncbi.gene.info$GeneID[match(ensembl$ensembl_gene_id[idx], ncbi.gene.info$ensembl_gene_id)]


####### HGNC

idx<- which(is.na(hgnc$ensembl_gene_id))
hgnc$ensembl_gene_id[idx] <- ensembl$ensembl_gene_id[match(hgnc$hgnc_id[idx], ensembl$hgnc_id)]

idx<- which(is.na(hgnc$entrez_id))
hgnc$entrez_id[idx] <-ncbi.gene.info$GeneID[match(hgnc$hgnc_id[idx], ncbi.gene.info$hgnc_id)]


####### NCBI

idx<- which(is.na(ncbi.gene.info$ensembl_gene_id))
ncbi.gene.info$ensembl_gene_id[idx] <- ensembl$ensembl_gene_id[match(ncbi.gene.info$GeneID[idx], ensembl$entrezgene_id)]

idx<- which(is.na(ncbi.gene.info$hgnc_id))
ncbi.gene.info$hgnc_id[idx] <- hgnc$hgnc_id[match(ncbi.gene.info$GeneID[idx], hgnc$entrez_id)]


sum(!is.na(ensembl$hgnc_id)) # 38616
sum(!is.na(ensembl$entrezgene_id)) # 26479

sum(!is.na(hgnc$ensembl_gene_id)) # 38420
sum(!is.na(hgnc$entrez_id)) # 41674

sum(!is.na(ncbi.gene.info$ensembl_gene_id)) # 26411
sum(!is.na(ncbi.gene.info$hgnc_id)) # 41674


gene1 <- hgnc$ensembl_gene_id[which(!is.na(hgnc$ensembl_gene_id))]
gene1

gene2 <- ensembl$ensembl_gene_id[which(!is.na(ensembl$hgnc_id))]
gene2

gene2[which(!gene2 %in% gene1)]


#############################################################################################
# ================================== Merge Annotations ==================================== #
#############################################################################################

####### HGNC + NCBI

hgnc.tmp <- hgnc[match(ncbi.gene.info$GeneID, hgnc$entrez_id),]
dim(hgnc.tmp)
dim(ncbi.gene.info)

idx <- which(!hgnc$hgnc_id %in% hgnc.tmp$hgnc_id)
idx

hgnc.tmp <- rbind(hgnc.tmp, hgnc[idx,])
dim(hgnc.tmp)
View(hgnc.tmp)

ncbi.gene.info[(nrow(ncbi.gene.info)+1):nrow(hgnc.tmp),] <- NA
dim(ncbi.gene.info)

ncbi.hgnc <- data.frame(ncbi.gene.info, hgnc.tmp)


which(ncbi.hgnc$GeneID!=ncbi.hgnc$entrez_id)
which(ncbi.hgnc$hgnc_id!=ncbi.hgnc$hgnc_id.1)
which(ncbi.hgnc$ensembl_gene_id!=ncbi.hgnc$ensembl_gene_id.1)

ncbi.hgnc[which(ncbi.hgnc$hgnc_id!=ncbi.hgnc$hgnc_id.1),]

idx <- which(is.na(ncbi.hgnc$ensembl_gene_id))
ncbi.hgnc$ensembl_gene_id[idx] <- ncbi.hgnc$ensembl_gene_id.1[idx]


####### HGNC + NCBI + ENSEMBL

ensembl.tmp <- ensembl[match(ncbi.hgnc$ensembl_gene_id, ensembl$ensembl_gene_id),]
dim(ensembl.tmp)

idx <- which(!ensembl$ensembl_gene_id %in% ensembl.tmp$ensembl_gene_id)
idx

ensembl.tmp <- rbind(ensembl.tmp, ensembl[idx,])
dim(ensembl.tmp)

ncbi.hgnc[(nrow(ncbi.hgnc)+1):nrow(ensembl.tmp),] <- NA
dim(ncbi.hgnc)

ncbi.hgnc.ensembl <- data.frame(ncbi.hgnc, ensembl.tmp)
View(ncbi.hgnc.ensembl)




###
idx <- which(is.na(ncbi.hgnc.ensembl$ensembl_gene_id))
idx
ncbi.hgnc.ensembl$ensembl_gene_id[idx] <- ncbi.hgnc.ensembl$ensembl_gene_id.2[idx]


###
idx <- which(is.na(ncbi.hgnc.ensembl$GeneID))
idx
ncbi.hgnc.ensembl$GeneID[idx] <- ncbi.hgnc.ensembl$entrez_id[idx]

idx <- which(is.na(ncbi.hgnc.ensembl$GeneID))
idx
ncbi.hgnc.ensembl$GeneID[idx] <- ncbi.hgnc.ensembl$entrezgene_id[idx]

###
idx <- which(is.na(ncbi.hgnc.ensembl$hgnc_id))
idx
ncbi.hgnc.ensembl$hgnc_id[idx] <- ncbi.hgnc.ensembl$hgnc_id.1[idx]

idx <- which(is.na(ncbi.hgnc.ensembl$hgnc_id))
idx
ncbi.hgnc.ensembl$hgnc_id[idx] <- ncbi.hgnc.ensembl$hgnc_id.2[idx]



#############################################################################################
# ===================================== Final Check ======================================= #
#############################################################################################

# ENSEMBL ID
test <- data.frame(ncbi.hgnc.ensembl$ensembl_gene_id, ncbi.hgnc.ensembl$ensembl_gene_id.1, 
                   ncbi.hgnc.ensembl$ensembl_gene_id.2, stringsAsFactors = F)
which(test$ncbi.hgnc.ensembl.ensembl_gene_id!=test$ncbi.hgnc.ensembl.ensembl_gene_id.1)
which(test$ncbi.hgnc.ensembl.ensembl_gene_id!=test$ncbi.hgnc.ensembl.ensembl_gene_id.2)

idx <- which(test$ncbi.hgnc.ensembl.ensembl_gene_id!=test$ncbi.hgnc.ensembl.ensembl_gene_id.1)
test[idx,]






# HGNC ID
test <- data.frame(ncbi.hgnc.ensembl$hgnc_id, ncbi.hgnc.ensembl$hgnc_id.1, 
                   ncbi.hgnc.ensembl$hgnc_id.2, stringsAsFactors = F)
which(test$ncbi.hgnc.ensembl.hgnc_id!=test$ncbi.hgnc.ensembl.hgnc_id.1)
which(test$ncbi.hgnc.ensembl.hgnc_id!=test$ncbi.hgnc.ensembl.hgnc_id.2)


# ENTREZ ID
test <- data.frame(ncbi.hgnc.ensembl$GeneID, ncbi.hgnc.ensembl$entrez_id, ncbi.hgnc.ensembl$entrez_id)
which(test$ncbi.hgnc.ensembl.GeneID!=test$ncbi.hgnc.ensembl.entrez_id)
which(test$ncbi.hgnc.ensembl.GeneID!=test$ncbi.hgnc.ensembl.entrez_id.1)


View(test)
