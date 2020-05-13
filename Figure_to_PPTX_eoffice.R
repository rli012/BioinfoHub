setwd('E:\\LABDATA/Haplotype/')

library(eoffice)


cvSummary <- read.table('Simulation.inferrenceSummary.all.txt', header = T, stringsAsFactors = F, sep = '\t')
cvSummary


GR <- as.matrix(cvSummary[,1:8])
GR


GRDa <- data.frame(missing = rep(seq(0,70,10),each=nrow(GR)),
                   pol=rep(3:15,8*8), 
                   fail = as.numeric(GR),
                   methods=factor(rep(rep(c('Hapi','PHMM'), each=13),8*4), 
                                  levels=c('Hapi', 'PHMM')),
                   snp = factor(rep(rep(c('5K','10K','50K','100K'),each=13*2),8),
                                levels=c('5K','10K','50K','100K')),
                   labColor = ifelse(as.numeric(GR)==0, 'blue', 'black'))


library(ggplot2)

htmp <- ggplot(data=GRDa, aes(x=missing, y=pol, fill=fail)) + 
  geom_tile(color='black') + 
  scale_fill_gradient('Incorrect inference',high='red', low='white', limits=c(0, 100),  na.value="gray",
                      breaks=seq(0,100,20),labels=seq(0,100,20))+ 
  facet_grid(methods~snp)+ # switch= 'both'
  labs(x="Missing genotype rate (%)", y="Number of gametes") + theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank()) +
  theme(legend.position='right', legend.title = element_blank()) + theme(legend.key.size=unit(0.2, "cm")) +
  theme(legend.key.width=unit(0.5, "cm"), legend.key.height = unit(1,'cm')) +
  scale_y_reverse(limits = rev(levels(GRDa$pol)), breaks=seq(3,15,1)) + 
  scale_x_continuous(breaks=seq(0,70,10)) + #position= top
  geom_text(data=GRDa,aes(label = fail), size=4) +
  theme(axis.ticks = element_blank(),
        axis.text=element_text(size=14, colour = 'black'),
        axis.title=element_text(size=18, colour = 'black')) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=16))

htmp

topptx(htmp, filename = 'E:\\Publications/Heatmap.pptx', width=11, height = 7)

