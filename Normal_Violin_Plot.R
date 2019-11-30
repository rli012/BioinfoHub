
pValue <- degFinal[rownames(sixGenes),'adj.P.Val']


my_annotations <- as.character(symnum(pValue, #corr = FALSE, na = FALSE,
                                      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***",'**','*','ns')))
my_annotations <- paste('FDR =', formatC(pValue, format = 'e', digits = 2))
my_annotations


#my_annotations <- ifelse(pValue >= 0.01, paste0('p = ', formatC(pValue, digits = 2)),
# paste0('p = ', formatC(pValue, format = "e", digits = 2)))
#anno <- data.frame(x=1:5,y=17.5,label=my_annotations)
library(dplyr)

exprMax <- dataForBoxPlot %>% group_by(gene) %>% summarise(maxExpr=max(expr))

anno <- data.frame(x=1.5,y=10.25,label=my_annotations,gene=exprMax$gene)
anno <- anno[which(anno$label!='ns'),]

p <- ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_violin(aes(fill=group),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_boxplot(aes(color=group, fill=group),
  #             outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
  #             outlier.fill = NA, alpha=1, width=0.5) +
  #stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
  #             fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_wrap(~gene, nrow=1) +
  geom_jitter(size=1, width=0.2, color='black') +
  labs(x='', y=expression('Expression Level (Log'[2]*'CPM)')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  geom_text(data =anno, aes(x, y, label=label, group=NULL),
            size=4) +
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size=12))+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_blank(),
        axis.title = element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p

