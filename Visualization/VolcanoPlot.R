library(ggplot2)


##################################################################
#========================= VolcanoPlot ==========================#
##################################################################

dataForVolcanoPlot <- readRDS(file='data/dataForVolcanoPlot.rds')
dataForVolcanoPlot

logFcThreshold <- log(2)
adjPvalThreshold <- 0.05

dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, logFC < logFcThreshold | FDR > adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, logFC >= logFcThreshold & FDR <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot, logFC <= -logFcThreshold & FDR <= adjPvalThreshold)] <- 'DOWN'
dataForVolcanoPlot$Significance

ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(FDR))) +
  #xlim(-7.5,7.5) 
  xlab(expression('log'[2]*'(Fold Change)')) + ylab(expression('-log'[10]*'(FDR)')) +
  geom_point(aes(color=Significance), alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold), 
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), color='darkgreen',linetype='dashed')+
  scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) + 
  scale_color_manual(values = c('green3',"black", "red")) + 
  geom_text_repel(data = subset(dataForVolcanoPlot, FDR < 10^-40 & logFC > log(2)), segment.alpha = 0.4,
                  aes(label = symbol), size = 3.5, color='red', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot, FDR < 10^-40 & logFC < log(2)*-1), segment.alpha = 0.4,
                  aes(label = symbol), size = 3.5, color='green3', segment.color = 'black') +
  #geom_text_repel(data = subset(degForVolcanoPlot, p<0.05 & log2FC<0), segment.alpha = 0.4,
  #                aes(label = Symbol), size = 3.5, color='green3', segment.color = 'black') +
  #geom_text_repel(data = subset(degForVolcanoPlot, (p<0.05 & log2FC>log2(1.2)) | (p<0.01 & log2FC>0)), segment.alpha = 0.4,
  #                aes(label = Symbol), size = 3, color='red', segment.color = 'black') +
  #geom_text_repel(data = subset(degForVolcanoPlot, (p<0.05 & log2FC<log2(1.2)*-1) | (p<0.01 & log2FC<0)), segment.alpha = 0.4,
  #                aes(label = Symbol), size = 3, color='green3', segment.color = 'black') +
  #geom_text_repel(data=subset(pairwiseComparisonStatistics.forPlot, pairwiseComparisonStatistics.forPlot$Significance %in% c('UP','DOWN')), 
  #                aes(label=Gene), segment.alpha = 0.4,size = 3, color='blue') +
  #geom_text_repel(data = subset(comparisonData, abs(ddCt) >= log2threshold & -log10(Pvalue) >= -log10(pvalThreshold)), segment.alpha = 0.4,
#                aes(label = Gene), size = annotSize) +
theme_gray() +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16),
        strip.text = element_text(size=14)) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
#theme_bw() +
#theme(axis.line = element_line(colour = "black"),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.border = element_rect(colour='black'),
#    panel.background = element_blank()) #+
#ggtitle(dataset)
#theme(axis.text=element_text(size=14), 
#      axis.title=element_text(size=16),
#      plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
#      legend.text = element_text(size = 14),
#      legend.title = element_text(size = 14, face = 'bold'),
#      strip.text = element_text(size = 14, face = 'bold'))
#scale_y_continuous(expand = c(0.3, 0)) +
#facet_wrap(~Comparison, ncol = 2) +
