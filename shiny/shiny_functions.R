google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'


dataset <- read_xlsx('data/PCa_Datasets.xlsx', sheet = 'shiny')

gene.annotation <- readRDS('data/Annotation/Homo_Sapiens_Gene_Annotation_ENSEMBL_HGNC_ENTREZ.RDS')
gene.annotation$alias_symbol <- gsub('"', '', gene.annotation$alias_symbol, fixed=T)

gene.default <- 'ENSG00000142515' # KLK3


boxplotFun <- function(dataForBoxPlot) {
  p <- ggplot(dataForBoxPlot, aes(x=group, y=expr)) + 
    geom_violin(aes(fill=group)) +
    #geom_boxplot(aes(fill=sample), width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #geom_boxplot(fill='white', width=0.2,
    #             outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
    #             outlier.fill = NA) +
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='white', fill='white') +
    #facet_wrap(~cell, nrow=1) +
    #geom_jitter(size=0.1, width=0.2) +
    #ylim(-0.1,3)+
    xlab('') + ylab(expression('log'[2]*'CPM')) +
    #ggtitle(paste0('Expression of ', gene.symbol)) +
    guides(fill = guide_legend(nrow=1)) +
    theme_bw() +
    theme(legend.position = 'bottom') +
    #theme(plot.title = element_text(hjust = 0.5, face='bold', size=16)) +
    theme(axis.text.y = element_text(size=12,color='black'),
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.spacing.x = unit(0.1, "cm"),
          axis.title = element_text(size=14),
          strip.text = element_text(angle=90, size=12, face='bold'),
          panel.border = element_rect(colour = "black"))
  
  p <- p + geom_jitter(size=0.1, width=0.2)
  
  return(p)
  
}

kmplotFun <- function(meta.data, expr.data, group, gene.id, quan) {
  
  if (group == 'time_to_death') {
    status <- 'os_status'
    x.title <- 'Overall Survival (months)'
  } else if (group == 'time_to_bcr') {
    status <- 'bcr_status'
    x.title <- 'Relapse-free Survival (months)'
  }  else if (group == 'time_to_metastasis') {
    status <- 'metastasis_status'
    x.title <- 'Metastasis-free Survival (months)'
  } 
  
  keep <- which(meta.data[,'sample_type'] %in% c('Primary','Tumor'))
  
  expr <- expr.data[gene.id,keep]
  
  daysToDeath <- meta.data[keep, group]
  vitalStatus <- meta.data[keep, status]
  
  expr.thres <- quantile(expr, quan/100, na.rm = T) # median(expr,na.rm=T) ##############
  expr.group <- expr > expr.thres
  
  dataForKMPlot <- data.frame(daysToDeath,vitalStatus, expr.group)
  
  nH <- sum(expr.group)
  nL <- sum(!expr.group)
  
  sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ expr.group)
  pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                          lower.tail = FALSE),digits=4)
  
  pValue <- format(1-pchisq(sdf$chisq, df=1),digits=3)
  
  HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  HR <- format(HR, digits = 3)
  upper95 <- format(upper95, digits = 3)
  lower95 <- format(lower95, digits = 3)
  
  label1 <- paste('HR = ', HR, ' (', lower95, '-', upper95, ')', sep='')
  label2 <- paste('P value = ', pValue, sep='')
  
  fit <- survfit(Surv(daysToDeath, vitalStatus) ~ expr.group, data=dataForKMPlot)
  
  lgdXpos <- 1/1.4
  lgdYpos = 0.9
  
  xpos = max(daysToDeath, na.rm=TRUE)/2.2
  ypos1 = 0.95
  
  p <- ggsurvplot(fit, data=dataForKMPlot, pval = paste(label1, '\n', label2), pval.coord = c(xpos, ypos1),
                  pval.size=4,
                  font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
                  #title = project,
                  #legend = c(lgdXpos, lgdYpos), 
                  #color = c('blue', 'green'),
                  palette= c('blue', 'red'),
                  legend.labs = c(paste('Low Expr (N=',nL,')',sep=''), 
                                  paste('High Expr  (N=',nH,')',sep='')),  
                  legend.title='group',
                  xlab = x.title, ylab = 'Survival probability',
                  #xlab = paste(type,'(months)'), ylab = 'Survival probability',
                  font.x = c(16), font.y = c(16), ylim=c(0,1), #16
                  ggtheme = theme_bw()+ theme(axis.line = element_line(colour = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              #panel.border = element_rect(colour='black'),
                                              panel.border = element_blank(),
                                              panel.background = element_blank(),
                                              legend.text = element_text(size=12),#14
                                              legend.title = element_blank(),
                                              legend.position = 'top',
                                              axis.text = element_text(size=14, color='black'))) #+
  
  return(p)
  
}


volcanoPlotFun <- function(dataForVolcanoPlot, logFcThreshold, adjPvalThreshold) {
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
      #geom_text_repel(data = subset(dataForVolcanoPlot, 
      #                              adj.P.Val < adjPvalThreshold & logFC > logFcThreshold), 
      #                segment.alpha = 0.4, aes(label = Symbol), 
      #                size = 3.5, color='red', segment.color = 'black') +
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
      
      return (p)
}



corrplotFun <- function(dataForCorrPlot) {
  
  dataForCorrPlot$group <- log2(dataForCorrPlot$group+1)
  
  corr <- cor.test(dataForCorrPlot$expr, dataForCorrPlot$group)$estimate
  p <- cor.test(dataForCorrPlot$expr, dataForCorrPlot$group)$p.value

  anno_text <- data.frame(
    label = paste0('corr = ', round(corr,3), '\n',
                   'p = ', ifelse(p >= 0.01,
                                  formatC(p, digits = 2),
                                  formatC(p, format = "e", digits = 2))),
    x     = max(dataForCorrPlot$expr),
    y     = max(dataForCorrPlot$group)
  )
  
  p <- ggplot(data=dataForCorrPlot, aes(x=expr, y=group)) +
    geom_point(size=1, color='darkred') +
    geom_smooth(method='lm') +
    #geom_text(data    = anno_text,
    #          mapping = aes(x = x, y = y, label = label),
    #          size=4.5) +
    #facet_wrap(~platform) +
    labs(x='Expression Level', y='Preoperative PSA') +
    theme_bw() +
    theme(axis.line = element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = 'none')+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=16),
          strip.text.x = element_text(size=14, face='bold'))
  
  return (p)
  
}


pieplotFun <- function(dataForPiePlot) {
  
  p <- ggplot(dataForPiePlot, aes(x = "", y = num, fill = sam)) +
    geom_bar(width = 1, stat = "identity", color = "white", size=0.5) +
    scale_fill_manual(values = c(google.yellow, google.blue, google.red, google.green)) +
    coord_polar("y", start = 0) + 
    geom_text(aes(label = num), position = position_stack(vjust = 0.5), size=4.5) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size=12)) +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
  
  return (p)
}

barplotFun <- function(dataForBarPlot) {
  
  p <- ggplot(data=dataForBarPlot, mapping=aes(x=group, y=num, fill=google.blue)) +
    geom_bar(stat='identity') +
    scale_x_discrete(limits=dataForBarPlot$group) +
    labs(x='', y='Number of samples') + 
    scale_fill_manual(values=c(google.blue)) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=12, color='black'),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title=element_text(size=14)) +
    theme(legend.text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = 'none') +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
  
  return (p)
}


histogramFun <- function(dataForHistogram) {

  p <- ggplot(data=dataForHistogram, aes(log2(preop_psa+1), fill=google.blue)) + 
    geom_histogram() + 
    scale_fill_manual(values=c(google.blue)) +
    labs(x=expression('Log'[2]*'(Preoperative PSA + 1)'), y='Count') +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     #panel.grid.major = element_blank(),
                     #panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
    theme(axis.text=element_text(size=12, color='black'),
          axis.title=element_text(size=14)) +
    theme(legend.text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = 'none') +
    theme(plot.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"))
    
  return (p)
  
  
}


