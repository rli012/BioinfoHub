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
                       
                       
### ComplexHeatmap
# have a symmetric continuous legend (say, from -3 to 3) with white in the center (at 0)
# https://github.com/jokergoo/ComplexHeatmap/issues/82
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# margin
# https://github.com/jokergoo/ComplexHeatmap/issues/174
column_names_max_height = unit(10, 'cm')
draw(ht, padding = unit(c(2, 2, 2, 20), "mm")) #bottom, left, top, right paddings
                       
                       
                       
### scatter
p1 <- ggplot(dataForScatterPlot, aes(x=x, y=y)) + geom_point(aes(color=Reg), size=1.5) + ggtitle('') +
  #geom_point(data = subset(dataForScatterPlot, Reg == 'UP'),
  #           aes(x, y, color=Reg), size=0.5) +
  #geom_point(data = subset(dataForScatterPlot, Reg == 'DOWN'),
  #           aes(x, y, color=Reg), size=0.5) +
  scale_color_manual(values = c('blue',"gray", "red")) +
  geom_abline(intercept = 1, slope=1, color='black', linetype='dashed') +
  geom_abline(intercept = -1, slope=1, color='black', linetype='dashed') +
  labs(x='x',y='y') +
  
  #geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label1.top), 
  #                aes(label=gene), segment.alpha = 0.4,size = 3, color='red',
  #                nudge_y = 1) +
  #geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label2.top), 
  #                aes(label=gene), segment.alpha = 0.4,size = 3, color='darkgreen',
  #                nudge_x = 1) +
  theme(legend.position = 'none') +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'none') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=14)) +
  theme(plot.title = element_text(hjust = 0.5, size=16, face='bold')) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) #+
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))

p1

                       
####### pheatmap
library(pheatmap)
library(colorspace)
col_fun = diverge_hcl(100, c = 100, l = c(50,90), power = 1.5)
                   
pheatmap(dataForHeatmap,
         scale = 'none',
         cluster_cols = F, 
         display_numbers = annoMatrix,
         border_color = NA,
         cluster_rows = F,
         fontsize_col = 12,
         fontsize_number = 10,
         treeheight_row = 0,
         show_rownames = T,
         annotation_legend = F,
         #breaks = c(seq(-mx,mx, mx*2/100)),
         breaks = c(seq(-8,8, 8*2/100)),
         color=col_fun
)

### Correlation

xpos <- (min(dataForCorrPlot$signature)+max(dataForCorrPlot$signature))/2
ypos <- as.numeric(summary(dataForCorrPlot$virus)[6])+0.5

p <- 0.049
c <- 0.458

ggplot(dataForCorrPlot, aes(x=signature, y=virus)) + 
  geom_point(aes(shape=group, color=group)) + 
  labs(x='', y='')+
  geom_smooth(method="lm",se=FALSE, col='darkgreen', size=0.5) + 
  scale_colour_manual(breaks = group, 
                      values = c('chocolate1', 'blue')) +
  ggplot2::annotate("text", x = xpos, y = ypos, 
                    label = paste0('Pearson Correlation = ', c, '\nP Value = ', p), size = 4.5) +
  theme_bw()+theme(legend.title = element_blank(),
                   legend.text = element_text(size=14),
                   axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='white'),
                   panel.background = element_blank(),
                   axis.text = element_text(size=14),
                   axis.title = element_text(size=16))

                       
### heatmap
library(ComplexHeatmap)
ht <- Heatmap(as.matrix(dataForHeatmap),
              #name = 'Expression',
              
              # COLOR
              #col = colorRampPalette(rev(c("red",'white','blue')), space = "Lab")(100),
              col=col_fun,
              na_col = 'grey',
              
              # MAIN PANEL
              column_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_dend = TRUE,
              column_order = o,
              show_column_dend = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE, 
              row_names_side = 'left',
              column_names_rot = 90,
              #column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 10),
              row_names_max_width = unit(15,'cm'),
              column_split = factor(c(rep("Module 1",8), rep('Module 2', 8), rep('Module 3', 9),
                                      rep('Module 4', 10), rep('Module 5', 8))),

              #column_names_max_height = unit(10, 'cm'),
              #column_split = factor(phenoData$Day,
              #                      levels=str_sort(unique(phenoData$Day), numeric = T)),
              
              #column_order = rownames(phenoData),
              
              # ANNOTATION
              #top_annotation = topAnnotation,
              
              
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(annoMatrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black", fill = 'black'))
              },
              
              # LEGEND
              heatmap_legend_param = list(
                at = c(-1, -0.5, 0, 0.5, 1),
                #labels = c("low", "zero", "high"),
                title = expression('Pearson Correlation'),
                title_position = 'leftcenter-rot',
                legend_height = unit(3, "cm"),
                just = c("right", "top")
              ))

draw(ht,annotation_legend_side = "right",row_dend_side = "left", heatmap_legend_side = "right")


                      
                       
####
                       
ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_boxplot(aes(fill=NULL),
               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
               outlier.fill = NA) +
  #geom_boxplot(aes(color=group, fill=group),
  #             outlier.shape = 21, outlier.size = 1,#outlier.colour = 'black',
  #             outlier.fill = NA, alpha=1, width=0.5) +
  #stat_summary(geom = "crossbar", width=0.45, fatten=0, color="white", position=p,
  #             fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  #ylim(-0.5,0.5)+
  facet_wrap(~signature, nrow=1) +
  geom_jitter(aes(color=group), size=2, width=0.05) +
  scale_color_manual(values = c(google.blue, google.red, google.green)) +
  labs(x='', y=expression('Fibroid Signature Score')) +
  #geom_segment(data=df,aes(x = x1, y = y1, xend = x2, yend = y2)) +
  #geom_text(data =anno, aes(x, y, label=label, group=NULL),
  #          size=4) +
  geom_signif(annotations = my_annotations, # optional
              comparisons = my_comparisons,
              step_increase = 0.1,
              vjust=.2,
              colour='gray20',
              tip_length=0.015) + # 0
  theme_bw()+
  theme(legend.position = 'none', legend.title = element_blank())+
  theme(axis.text = element_text(size=12,color='black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size=16),
        strip.text = element_text(size=14, face='bold')) +
  theme(panel.border = element_rect(color='grey60')) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
