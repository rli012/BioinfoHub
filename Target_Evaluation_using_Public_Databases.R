library(ggplot2)
library(dplyr)


##### Gene Expression
### GTEx
### TCGA
### Human Protein Atlas
### Human Proteome Map

#=== For Immune Cells ===#
### DICE
### Blueprint


##### Genetic Association
### GWAS catalog
### Phenoscanner
### PheWAS
### ClinVar
### GWAS Central
### UK10K
### DisGeNET

### TOPMed
### ExAC



##### Knockout Phenotype
### MGI

##### Target-Disease Association
### DisGeNET
### Open Targets



exprData <- read.table('data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', header = T, stringsAsFactors = F, 
                       sep = '\t', comment.char = '#', skip = 2)
head(exprData)

exprData[1:5,1:5]

genes <- sapply(exprData$Name, function(x) strsplit(x, '.', fixed = T)[[1]][1])
genes

filter <- which(duplicated(genes))
filter

exprData <- exprData[-filter,]

rownames(exprData) <- sapply(exprData$Name, function(x) strsplit(x, '.', fixed = T)[[1]][1])



phenoSubject <- read.table('data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', header = T, stringsAsFactors = F, sep = '\t')
phenoSubject

phenoSample <- read.table('data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', header = T, stringsAsFactors = F, sep = '\t',
                          quote = '')
phenoSample


phenoSample$PATID <- paste(sapply(phenoSample$SAMPID, function(x) strsplit(x, '-')[[1]][1]),
                           sapply(phenoSample$SAMPID, function(x) strsplit(x, '-')[[1]][2]),
                           sep='-')

idx <- match(phenoSample$PATID, phenoSubject$SUBJID)
idx


phenoSample$SEX <- phenoSubject$SEX[idx]
phenoSample$AGE <- phenoSubject$AGE[idx]
phenoSample$DTHHRDY <- phenoSubject$DTHHRDY[idx]

phenoSample$SEX[phenoSample$SEX==1] <- 'Male'
phenoSample$SEX[phenoSample$SEX==2] <- 'Female'

phenoData <- phenoSample


rownames(phenoData) <- gsub('-', '.', phenoData$SAMPID, fixed = T)

ovlp <- intersect(rownames(phenoData),colnames(exprData))
length(ovlp)


phenoData <- phenoData[ovlp,]
exprData <- exprData[,ovlp]

rownames(exprData)

saveRDS(exprData, 'data/GTEx_exprData.RDS')
saveRDS(phenoData, 'data/GTEx_phenoData.RDS')

###########################################################################################

g <- 'ENSG0000010xxxx'
#idx <- which(rownames(exprData)==g)
idx
expr <- exprData[g,]
expr


dataForBoxPlot <- data.frame(expr=log2(as.numeric(expr)+0.1), tissue=phenoData$SMTS, subtype=phenoData$SMTSD, gender=phenoData$SEX, age=phenoData$AGE)
dataForBoxPlot

o <- dataForBoxPlot %>% group_by(tissue) %>%
  summarise(med=median(expr))

o <- as.character(o$tissue[order(o$med, decreasing = F)])
o


dataForBoxPlot$tissue <- factor(dataForBoxPlot$tissue, levels = o)
dataForBoxPlot$tissue

ggplot(data=dataForBoxPlot, aes(x=tissue, y=expr)) +
  geom_boxplot(aes(fill=tissue),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_boxplot(aes(color=group, fill=group),
  #             outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
  #             outlier.fill = NA, alpha=1, width=0.5) +
  #stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
  #             fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  #facet_wrap(~tissue, nrow=1) +
  geom_jitter(size=0.5, width=0.2, color='darkblue') +
  labs(x='', y=expression('Expression Level (Log'[2]*'TPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=5) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


###
for (tissue in o) {
  print (tissue)
  expr <- dataForBoxPlot$expr[dataForBoxPlot$tissue==tissue]
  group <- dataForBoxPlot$gender[dataForBoxPlot$tissue==tissue]
  
  if (length(unique(group))==1) {
    next
  }
  
  t <- t.test(expr~group)
  p <- t$p.value
  
  print (p)
}
  
  



dataForBarPlot <- dataForBoxPlot %>% group_by(tissue) %>%
  summarise(me=mean(expr), sd=sd(expr))

o <- order(dataForBarPlot$expr, decreasing = F)

dataForBarPlot$Subtype <- factor(dataForBarPlot$Subtype, levels=dataForBarPlot$Subtype[o])


# colors
ggplot(data=dataForBarPlot, aes(x=tissue, y=me, fill=tissue, color=tissue)) +
  geom_bar(stat='identity', width=.6) + #coord_flip()
  geom_errorbar(aes(ymin=ifelse(me>0, me, me-sd),
                    ymax=ifelse(me>0, me+sd, me)),
                width=.5, size=0.5,
                position=position_dodge(.9)) +
  labs(x='',y=expression('Expression Level (Log'[2]*'TPM)')) +
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  #scale_fill_manual(values = rep('black',nrow(dataForBarPlot))) +
  #scale_color_manual(values = rep('black',nrow(dataForBarPlot))) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'none') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))




##########
expr <- exprData[g,]
expr

dataForBoxPlot <- data.frame(expr=log2(as.numeric(expr)+0.1), tissue=phenoData$SMTS, subtype=phenoData$SMTSD, gender=phenoData$SEX, age=phenoData$AGE)
dataForBoxPlot

keep <- which(dataForBoxPlot$tissue %in% c('Liver','Ovary','Testis','Lung','Adrenal Gland', 'Breast'))
dataForBoxPlot <- dataForBoxPlot[keep,]


o <- dataForBoxPlot %>% group_by(tissue) %>%
  summarise(med=median(expr))

o <- as.character(o$tissue[order(o$med, decreasing = F)])
o


dataForBoxPlot$tissue <- factor(dataForBoxPlot$tissue, levels = o)
dataForBoxPlot$tissue


ggplot(data=dataForBoxPlot, aes(x=gender, y=expr)) +
  geom_boxplot(aes(fill=gender),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_boxplot(aes(color=group, fill=group),
  #             outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
  #             outlier.fill = NA, alpha=1, width=0.5) +
  #stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
  #             fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_wrap(~tissue, nrow=1) +
  geom_jitter(size=0.5, width=0.2, color='darkblue') +
  labs(x='', y=expression('Expression Level (Log'[2]*'TPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=5) +
  theme(legend.position = 'right')+
  theme(axis.text = element_text(size=14,color='black'),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title = element_text(size=16),
        strip.text = element_text(size=14, face='bold', angle = 90)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


dataForBoxPlot




p <- position_dodge(0.58)

bxplt <- ggplot(dataForBoxPlot, aes(x=tissue, y=expr, group = interaction(gender, tissue)))
#fill=interaction(group,project))
#stat_boxplot(geom ='errorbar', width=0.5, position = position_dodge(width = 0.75))
#stat_summary(geom ='crossbar', width=0.5, color='white')

bxplt+geom_boxplot(#aes(color=group, fill=group),
                   outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
                   outlier.fill = NA, alpha=1, width=0.5,
                   position = p) +
  stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  #scale_fill_npg(palette = 'nrc', alpha = 0.5)+
  #scale_fill_manual(values=c('limegreen', 'blue'),
  #                  labels=c('Normal', 'Dysplasia'),
  #                  name='Sample type') +
  #geom_jitter(size=0.3)+
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #ylim(-5,20) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=5) +
  #geom_signif(annotations = my_annotations,
  #            y_position = rep(7.5,5),
  #            xmin = c(1,2,3,4,5)-0.15,
  #            xmax = c(1,2,3,4,5)+0.15,
  #            #comparisons = my_comparisons,
  #            step_increase = 0,
  #            vjust=.2,
  #            colour='gray20',
  #            tip_length=0) + # 0
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 0, hjust=0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))





ggplot(dataForBoxPlot, aes(x=tissue, y=expr, fill=gender)) + # fill=interaction(group,project))
  stat_boxplot(geom ='errorbar', width=0.5, position = position_dodge(width = 0.75))+
  
  geom_boxplot(outlier.colour = 'black', outlier.shape = 21,
               outlier.fill = NA) +
  scale_fill_manual(values=c('limegreen', 'blue'),
                    labels=c('Tumor', 'Normal'),
                    name='Sample type') +
  xlab('')+ylab('Expression value (log2CPM)') +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'top') +
  theme(axis.title=element_text(size=16), 
        axis.text = element_text(color='black', size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color='black'),
        panel.background = element_blank())


