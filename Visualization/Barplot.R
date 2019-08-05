

barData <- readRDS(file='data/dataForBarPlot_Errorbar.rds')
barData

dataForBarPlot <- data.frame(expr=apply(barData, 1, function(x) mean(x, na.rm=T)),
                             sd = apply(barData, 1, function(x) sd(x, na.rm=T)),
                             cell=rownames(barData))

o <- order(apply(barData, 1, function(x) mean(x, na.rm=T)), decreasing = T)

dataForBarPlot$cell <- factor(dataForBarPlot$cell, levels=dataForBarPlot$cell[o])


ggplot(data=dataForBarPlot, aes(x=cell, y=expr), fill='black', color='black') + 
  geom_bar(stat='identity', width=.8) + #coord_flip()
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y='mRNA Level (TPM)') + 
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
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.5, unit = "cm"))

