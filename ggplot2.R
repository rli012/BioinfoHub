### Diagnal, abline
ggplot(data=dataForScatterPlot, aes(x=, y)) +
  geom_point(size=2, color='darkred') +
  geom_abline(intercept = 0, slope=1, color='blue', linetype='dashed')+
  xlim(6.5,13.5) + ylim(7,13.5) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text.x = element_text(size=14, face='bold'))


### Density plot
ggplot(dataForDensityPlot, aes(x=expr, fill=array)) +
  geom_density(alpha=0.4) + #xlim(3,12) +
  #geom_vline(xintercept=c(10,20), linetype='dashed', size=0.5) +
  labs(x=expression('Log'[2]*'Intensity'), y='Density') +
  theme_gray() +
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=12)) +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16),
        strip.text = element_text(size=14))
