setwd('C:\\Users/rli3/Documents/miRNomes/')
library(limma)
library(ggplot2)
library(edgeR)
library(ggrepel)


google.red <- '#ea4235'
google.yellow <- '#fabd03'
google.green <- '#34a853'
google.blue <- '#4286f5'

meta.tcga <- readRDS('shinyApp/data/Metadata_TCGA.RDS')

mir.tcga <- readRDS('shinyApp/data/miRNA_Expression_TCGA.RDS')

mir.annotation <- readRDS('shinyApp/data/miRBase_10.0_22.RDS')


meta <- meta.tcga[['TCGA-PRAD']]
expr <- mir.tcga[['TCGA-PRAD']]



expr[1:5,1:5]

colnames(expr) == rownames(meta)


group <- meta$sample_type
group


####

filter <- which(meta$sample_type=='Normal')
filter

meta <- meta[-filter,]
expr <- expr[,-filter]

group <- sapply(meta$gleason_score, function(x) strsplit(x, ' ')[[1]][1])
group <- as.numeric(group)
group

filter <- which(is.na(group))
filter


group <- group[-filter]
group

table(group)
group <- ifelse(group>=8, 'High', 'Low')
group



group <- meta$pathologic_T
group

filter <- which(group=='NA')
filter

group <- group[-filter]

group <- ifelse(grepl('T2',group), 'Low', 'High')

deg.group <- factor(group)

design <- model.matrix(~0+deg.group)
colnames(design) <- levels(deg.group)

contrast.matrix <- makeContrasts(contrasts='High - Low',
                                 levels=design)
contrast.matrix

### Differential gene expression analysis (limma)

fit <- lmFit(expr[,-filter], design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')

dataForVolcanoPlot <- dgeTable

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.01

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC < logFcThreshold | adj.P.Val > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC >= logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, 
                                     logFC <= -logFcThreshold & adj.P.Val <= adjPvalThreshold)] <- 'DOWN'
dataForVolcanoPlot$miRNA.ID

dataForVolcanoPlot <- data.frame(miRNA.Accession=rownames(dataForVolcanoPlot),
                                 miRNA.ID=mir.annotation[rownames(dataForVolcanoPlot),]$Name,
                                 dataForVolcanoPlot,
                                 stringsAsFactors = F)

p <- ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(adj.P.Val))) +
  #xlim(-2,2) +
  labs(x=expression('log'[2]*'(Fold Change)'), 
       y=(expression('-log'[10]*'(FDR)')), 
       title=NULL) +
  geom_point(aes(color=Significance), alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  scale_color_manual(values = c(google.green,"gray", google.red)) +
  #facet_wrap(~Comparison, ncol = 2) +
  geom_point(data = subset(dataForVolcanoPlot,
                           miRNA.ID %in% c('hsa-miR-7-1-3p','hsa-miR-7-5p')),
                  aes(x=logFC, y=-log10(adj.P.Val)), size=2, color='blue', fill='gray',shape=21)+
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                miRNA.ID %in% c('hsa-miR-7-1-3p','hsa-miR-7-5p')),
                 aes(label = miRNA.ID),
                 segment.alpha = 1, segment.size = 0.5, 
                 min.segment.length = 5,
                 
                 size = 4, color='black', segment.color = 'black') +
  
  
  #geom_text_repel(data = subset(dataForVolcanoPlot, 
  #                              adj.P.Val < adjPvalThreshold & logFC < logFcThreshold*-1), 
  #                segment.alpha = 0.4, aes(label = Symbol), 
  #                size = 3.5, color='green3', segment.color = 'black') +
  
  theme_bw() +
  theme(axis.line = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))


p




group <- sapply(meta$gleason_score, function(x) strsplit(x, ' ')[[1]][1])
#group <- as.numeric(group)
group

filter <- which(group=='NA')
filter

group <- group[-filter]
group

group <- factor(group, levels=c('6','7','8','9','10'))



group <- meta$pathologic_T
group

filter <- which(group=='NA')
filter

group <- group[-filter]



dataForBoxPlot <- data.frame(expr=c(expr['MIMAT0004553',-filter], expr['MIMAT0000252',-filter]),
                             group=rep(group,2),
                             gene=rep(c('hsa-miR-7-1-3p','hsa-miR-7-5p'), each=length(group)),
                             stringsAsFactors = F)

                             
dataForBoxPlot$gene <- factor(dataForBoxPlot$gene,
                               levels=c('hsa-miR-7-1-3p','hsa-miR-7-5p'))


ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_boxplot(aes(fill=group),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_smooth(method='lm') +
  facet_wrap(~gene) +
  geom_jitter(size=2, width=0.05, color='black') +
  labs(x='Gleason Score', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  # geom_signif(annotations = my_annotations, # optional
  #             comparisons = my_comparisons,
  #             step_increase = 0.1,
  #             vjust=.2,
  #             colour='gray20',
  #             tip_length=0.015) + # 0
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16)) +
  theme(strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




#################################################################


meta <- readRDS('shinyApp/data/fromTCGA/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')
expr <- mir.tcga[['TCGA-PRAD']]

samples <- intersect(rownames(meta), colnames(expr))
samples

meta <- meta[samples,]
expr <- expr[,samples]


filter <- which(meta$TherapiesYesNo=='YES')
filter

meta <- meta[-filter,]
expr <- expr[,-filter]

table(meta$gleason_score)


colnames(expr) == rownames(meta)

meta$rfs <- ifelse(meta$recurrence_status==1, 
                   meta$days_to_first_biochemical_recurrence,
                   meta$days_to_last_followup)

meta$rfs

os.time <- meta$rfs/365*12
os.status <- meta$recurrence_status
os.status

keep <- which(os.time < 60)
keep

idx <- which(os.time>60)
idx

os.time[idx] <- 60
os.status[idx] <- 0

# os.time <- meta$OS.time/365*12
# os.status <- meta$OS
# 
# os.time



# os.time <- as.numeric(meta$days_to_first_biochemical_recurrence)/365*12
# os.time
# 
# nonComplt <- is.na(os.time)
# 
# os.status <- as.numeric(ifelse(nonComplt, 0, 1))
# os.time[nonComplt] <- as.numeric(meta$days_to_last_followup[nonComplt])/365*12
# os.time


mir.id <- 'MIMAT0004553'
mir.id <- 'MIMAT0000252'
mir.name <- mir.annotation[mir.id,]$Name

# os.time <- as.numeric(meta$OS.time)/30
# os.status <- as.numeric(meta$OS)

score <- expr[mir.id,]
coxtest <- coxph(Surv(os.time, os.status) ~ score)
summcph <- summary(coxtest)

coeffs <- c(summcph$coefficients[,2], summcph$conf.int[,3:4], 
            summcph$coefficients[,5])

coeffs


###
risk.group <- score < median(score, na.rm = T)

n.high <- sum(risk.group, na.rm=T)
n.low <- sum(!risk.group, na.rm=T)

sdf <- survdiff(Surv(os.time, os.status) ~ risk.group)
p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
#p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

coeffs <- c(hr, lower95, upper95, p.val)
coeffs


dataForKMPlot <- data.frame(expr=score, os.time=os.time, 
                            os.status=os.status)

# dataForKMPlot <- data.frame(expr=score[keep], os.time=os.time[keep], 
#                             os.status=os.status[keep])


dataForKMPlot

p <- KMPlotFun(dataForKMPlot, type='rfs', sep = '3rdQu')
p




KMPlotFun <- function(dataForKMPlot, sep='median', type='os') {
  
  if (sep=='1stQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[2])
  } else if (sep=='median') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[3])
  } else if (sep=='mean') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[4])
  } else if (sep=='3rdQu') {
    risk.threshold <- as.numeric(summary(dataForKMPlot$expr)[5])
  }
  
  dataForKMPlot$risk.group <- dataForKMPlot$expr > risk.threshold
  
  if (type == 'os') {
    x.title <- 'Overall Survival (months)'
  } else if (type == 'rfs') {
    x.title <- 'Relapse-free Survival (months)'
  }  else if (type == 'mfs') {
    x.title <- 'Metastasis-free Survival (months)'
  }
  
  n.high <- sum(dataForKMPlot$risk.group, na.rm=T)
  n.low <- sum(!dataForKMPlot$risk.group, na.rm=T)
  
  sdf <- survdiff(Surv(dataForKMPlot$os.time, dataForKMPlot$os.status) ~ dataForKMPlot$risk.group)
  p.val <- pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE)
  #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  hr = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(hr) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(hr) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  hr <- format(hr, digits = 2, nsmall=2)
  upper95 <- format(upper95, digits = 2, nsmall=2)
  lower95 <- format(lower95, digits = 2, nsmall=2)
  
  p.val <- ifelse(p.val >= 0.01, formatC(p.val, digits = 2), 
                  formatC(p.val, format = "e", digits = 2))
  
  label.hr <- paste('HR = ', hr, ' (', lower95, ' - ', upper95, ')', sep='')
  label.p <- paste('P Value = ', p.val, sep='')
  
  fit <- survfit(Surv(os.time, os.status) ~ risk.group, data=dataForKMPlot)
  
  lgd.xpos <- 0.3
  lgd.ypos = 0.22
  
  p.xpos = max(dataForKMPlot$os.time, na.rm=TRUE)/50
  p.ypos = 0.05
  
  #title <- 'PFR10YR'
  #type <- 'Relapse-free Survival'
  
  plt <- ggsurvplot(fit, data=dataForKMPlot, pval = paste0(label.hr, '\n', label.p), pval.coord = c(p.xpos, p.ypos),
                    pval.size=4.2,
                    font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                    #title = title,
                    legend = c(lgd.xpos, lgd.ypos), 
                    #color = c('blue', 'green'),
                    palette= c(google.blue, google.red),
                    legend.labs = c(paste('Low Expression (N=',n.low,')',sep=''), 
                                    paste('High Expression (N=',n.high,')',sep='')),  
                    legend.title='', # Group
                    xlab = x.title, ylab = 'Survival Probability',
                    font.x = c(16), font.y = c(16), ylim=c(0,1), #20
                    ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                #panel.border = element_rect(colour='black'),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.text = element_text(size=12),#16
                                                legend.title = element_blank(), # 16
                                                legend.box.background = element_blank(),
                                                axis.title = element_text(size = 16, face = 'bold'),
                                                axis.text = element_text(size=12, color='black', face = 'bold'))) # 18
  
  return(plt[[1]])
  
}






