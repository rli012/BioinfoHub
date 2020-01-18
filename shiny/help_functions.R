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

