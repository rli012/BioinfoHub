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

### World map
library(maps)

world.map <- map_data('world')
summary(world.map$lat)
summary(world.map$long)

df <- data.frame(long=c(119.4101,-122.2711),
                 lat=c(35.9957,37.5585),
                 level=c(1,-1))

ggplot(world.map, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill='lightyellow', color="gray80", alpha=0.5, size=0.1) +
  #coord_quickmap() +
  #labs(x='',y='') +
  #ylim(0,50) +
  geom_point(data=df, aes(long, lat, color=level), size=2) +
  scale_colour_gradientn(limits=c(-1.5,1.5),
                         colors= c("blue",'white',"red")) + 
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   #panel.grid.minor = element_blank(),
                   panel.grid.minor = element_line(color='gray', linetype = 'dashed', size=0.2),
                   panel.grid.major = element_line(color='gray', linetype = 'dashed', size=0.2),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  #ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
  #theme(axis.text=element_text(size=14, color='black'),
  #      axis.text.x =element_text(size=14, color='black', angle=45, hjust=1),
  #      axis.title=element_text(size=15)) +
  theme(plot.title = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.text = element_text(size = 10), #v.just=-1
        legend.title = element_text(size = 12),
        #legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = c(0.05,0.3)) +
  theme(strip.text = element_text(size = 14),
        legend.key.size = unit(0.8,'cm'))
