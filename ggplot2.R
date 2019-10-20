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

### Volcano plot
ggplot(dataForVolcanoPlot, aes(x = logFC, y = -log10(P.Value))) +
  #xlim(-7.5,7.5)
  labs(x=expression('log'[2]*'(Fold Change)'), 
       y=(expression('-log'[10]*'(FDR)')), 
       title=NULL) +
  geom_point(color='#90909090', fill='#d0d0d090', alpha=1, size=1, shape=21) +
  geom_point(data=subset(dataForVolcanoPlot, Gene %in% genes),
                         color='navy', fill='#0000FFa0', alpha=1, size=2, shape=21) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed') +
  #scale_x_continuous(breaks=c(-4,-2,0,2,4,6,8,10)) +
  #scale_y_continuous(expand = c(0.3, 0)) +
  #scale_color_manual(values = c('#4285F4',"gray", '#FBBC05')) +
  #scale_color_manual(values = c('green3',"black", "red")) +
  #facet_wrap(~Comparison, ncol = 2) +
  geom_text_repel(data=subset(dataForVolcanoPlot, Gene %in% genes),
                  segment.alpha = 0.4, aes(label = Gene), 
                  size = 3.5, color='red', segment.color = 'black') +
  #geom_text_repel(data = subset(dataForVolcanoPlot, 
  #                              adj.P.Val < adjPvalThreshold & logFC > log2(2)), 
  #                segment.alpha = 0.4, aes(label = Symbol), 
  #                size = 3.5, color='red', segment.color = 'black') +
  #geom_text_repel(data = subset(dataForVolcanoPlot, 
  #                              adj.P.Val < adjPvalThreshold & logFC < log2(2)*-1), 
  #                segment.alpha = 0.4, aes(label = Symbol), 
  #                size = 3.5, color='green3', segment.color = 'black') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank()) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
        legend.position = 'none',
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14, face = 'bold'),
        strip.text = element_text(size = 14, face = 'bold'))


### World map
library(maps)

world.map <- map_data('world')
summary(world.map$lat)
summary(world.map$long)

df <- data.frame(long=c(119.4101,-122.2711),
                 lat=c(35.9957,37.5585),
                 level=c(1,-1))

ggplot(world.map, aes(long, lat)) +
  geom_map(map=world.map, aes(map_id=region), fill='lightyellow', color="gray80", alpha=0.5, size=0.1) +
  #coord_quickmap() +
  #labs(x='',y='') +
  #xlim(-10,28) + ylim(36,63) +
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

### Boxplot, geom_jitter, set colors
color <- ifelse(sapply(1:10, function(i) exprData[genes[i],] > mControl[i]), 'royalblue', 'lightblue')
color <- ifelse(sapply(1:10, function(i) exprData[genes[i],] > mControl[i]), 'darkblue', 'lightblue')

dataForBoxPlot <- data.frame(expr=as.numeric(t(exprData[genes,])), 
                             colors=c(color),
                             group=phenoData$Status,
                             gene=rep(genes, each=ncol(exprData)),
                             threshold=rep(mControl, each=ncol(exprData)),
                             stringsAsFactors = F)

dataForBoxPlot$gene <- factor(dataForBoxPlot$gene, levels=genes)
dataForBoxPlot$gene

threshold <- data.frame(gene=genes, thres=mControl, stringsAsFactors = F)
threshold$gene <- factor(threshold$gene, levels=genes)

ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_boxplot(aes(fill=NULL),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  geom_hline(data=threshold, aes(yintercept = thres), color='blue',linetype='dashed', size=0.5)+
  #ylim(6,9)+
  #geom_boxplot(aes(color=group, fill=group),
  #             outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
  #             outlier.fill = NA, alpha=1, width=0.5) +
  #stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
  #             fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_wrap(~gene, nrow=2) +
  geom_jitter(size=1.2, width=0.2, aes(colour=colors)) + 
  #scale_color_manual(values = dataForBoxPlot$colors) + # Jitter color palette
  labs(x=NULL, y=expression('Log'[2]*'Intensity')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=5) +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size=14,color='black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
                       
                       
#### heatmap for single cell
                       
levels(cell_type) <- seq(nlevels(cell_type)) #make cluster names short as numeric

col.Cluster <- list(`Cell type`=colorspace::qualitative_hcl(nlevels(cell_type), palette="Dark 3"))
names(col.Cluster$`Cell type`) <- levels(cell_type)
#col.Signatures <- list(Signatures=colorspace::qualitative_hcl(nlevels(sig.ord$CellType), palette="Dark 3"))
#names(col.Signatures$Signatures) <- levels(sig.ord$CellType)
col.Cluster

#annoColors <- list(
#  `Sample Type`=c(Normal='darkolivegreen',
#                  Dysplasia='lightcoral'))


#topAnnotation = HeatmapAnnotation(`Sample Type`=phenoData[,'SampleType'],
#                                  `Histology Grade`=phenoData[,'HistologyGrade'],
#                                  col=annoColors,
#                                  simple_anno_size_adjust = TRUE,
#                                  #annotation_height = c(1,1),
#                                  height = unit(8, "mm"),
#                                  #summary = anno_summary(height = unit(4, "cm")),
#                                  show_legend = c("bar" = TRUE),
#                                  show_annotation_name = F)


ht <- Heatmap(as.matrix(sig.exprs), name='Log2(UMI)', show_row_names = TRUE, show_column_names = FALSE, 
        top_annotation = HeatmapAnnotation(`Cell type` = cell_type, 
                                           annotation_legend_param = 
                                             list(`Cell type` = list(at = levels(cell_type),
                                                                 labels = levels(s[['RNA']]@misc$cell_type))),
                                           col=col.Cluster,
                                           show_legend=TRUE, 
                                           show_annotation_name=FALSE),
        #top_annotation = HeatmapAnnotation(Cluster = cell_type, col=col.Cluster, show_legend=FALSE, show_annotation_name=FALSE),
        #left_annotation = rowAnnotation(Signatures = sig.ord$CellType, col=col.Signatures, show_legend=FALSE, show_annotation_name=FALSE),
        #left_annotation = rowAnnotation(Signatures = sig.ord$CellType, show_legend=FALSE, show_annotation_name=FALSE),
        row_names_side = "left", row_names_gp = gpar(fontsize = 7), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = 'bold'),
        column_title_gp = gpar(fontsize = 8),
        #row_split = sig.ord$CellType,
        column_split = cell_type,
        cluster_rows = TRUE, cluster_columns = FALSE, 
        show_column_dend = FALSE, show_row_dend = FALSE,
        cluster_column_slices = FALSE,
        show_heatmap_legend = TRUE,
        #heatmap_legend_param = list(legend_direction='horizontal'),
        col = colorRamp2(c(0, floor(max(sig.exprs))), c("grey95", "red")))
ht

                       
### Pie
ggplot(dataForPiePlot, aes(x = "", y = Prop, fill = Tissue)) +
  geom_bar(width = 1, stat = "identity", color = "white", size=0.1) +
  coord_polar("y", start = 0) + 
  #geom_text(aes(y = lab.ypos, label = Prop), color = "white") +
  #geom_text(aes(label = paste0(Prop, "%")), position = position_stack(vjust = 0.5)) +
  #scale_fill_manual(values = mycols) +
  theme_void() + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
plot.margin = margin(t=0.1,l=0.1,r=0.1,b=0.1, unit = 'cm'))
