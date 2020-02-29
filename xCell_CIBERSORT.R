dge <-  DGEList(counts = countMatrix)
rpkmData <- rpkm(dge, gene.length = geneLength)

keep <- rowSums(edgeR::cpm(dge) > 1) >= 8 #0.5*ncol(countMatrix)
sum(keep)

rpkmData <- rpkmData[keep,]

write.table(rpkmData, file='report/RPKM_for_CYBERSORT.txt', quote=F, sep='\t')
write.table(rpkmData, file='report/RPKM_for_ssGSEA.txt', quote=F, sep='\t')



cibersort <- t(xCell[,rownames(phenoData)])
cibersort

idx <- which(phenoData$Visit %in% c('BASELINE'))
idx

dataForBoxPlot <- data.frame(expr=as.numeric(unlist(cibersort[idx,])),
                             group=rep(phenoData$sAg_loss[idx],ncol(cibersort)),
                             gene=rep(colnames(cibersort),each=nrow(phenoData[idx,])),
                             stringsAsFactors = F)
dataForBoxPlot



m <- dataForBoxPlot %>% group_by(gene, group) %>%
  summarise(mean(expr))

m

o <- order(m$`mean(expr)`[m$group=='sAg_loss'], decreasing = T)
o

genes <- m$gene[m$group=='sAg_loss'][o]
genes


s1 <- round(m$`mean(expr)`[m$group=='sAg_loss'][o],4)
s2 <- round(m$`mean(expr)`[m$group=='No_sAg_loss'][o],4)

s <- paste(s1,s2, sep='/')

p <- round(apply(cibersort[idx,], 2, function(v) wilcox.test(v~phenoData[idx,]$sAg_loss)$p.value),3)
p <- p[genes]

sort(p)

p <- paste0(s, ' (p=', p, ')')
p

names(p) <- genes

rowAnno = rowAnnotation(sig = anno_mark(at = 1:67, labels = p[rownames(t(cibersort[idx,]))]))
rowAnno


f <- m$gene[m$group=='sAg_loss'][o[1]]
f

o <- order(phenoData[idx,]$sAg_loss, cibersort[idx, f], decreasing = F)
o

samples <- rownames(phenoData[idx,])[o]
samples


annoColors <- list(
  `sAg_loss`=c(sAg_loss='darkolivegreen',
               No_sAg_loss='lightcoral'),
  `Visit`=c(BASELINE='gray96',
            WEEK12='gray48',
            WEEK24='gray0'))

topAnnotation = HeatmapAnnotation(`sAg_loss`=phenoData[idx,'sAg_loss'],
                                  #`Visit`=phenoData[idx,'Visit'],
                                  col=annoColors,
                                  simple_anno_size_adjust = TRUE,
                                  #annotation_height = c(1,1),
                                  height = unit(8, "mm"),
                                  #summary = anno_summary(height = unit(4, "cm")),
                                  show_legend = c("bar" = TRUE),
                                  show_annotation_name = F)


col_fun = colorRampPalette(rev(c("red",'yellow','blue')), space = "Lab")(100)
#col_fun = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
# col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")) # For Legend()
#library(colorspace)
#col_fun = diverge_hcl(100, c = 100, l = c(50,90), power = 1.5)


### pre-defined order of columns
ht <- Heatmap(as.matrix(t(cibersort[idx,])),
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
              show_column_dend = FALSE,
              show_row_names = TRUE,
              show_column_names = FALSE,
              column_names_rot = 90,
              column_names_gp = gpar(fontsize = 10),
              #column_names_max_height = unit(3, 'cm'),
              column_split = factor(phenoData[idx,]$sAg_loss),
              #levels=str_sort(unique(phenoData$Day), numeric = T)),
              
              column_order = samples,
              row_order = genes,
              
              row_names_side = "left",
              
              # ANNOTATION
              top_annotation = topAnnotation,
              
              right_annotation = rowAnno,
              
              # LEGEND
              heatmap_legend_param = list(
                #at = c(-5, 0, 5),
                #labels = c("low", "zero", "high"),
                title = "Score",
                title_position = 'leftcenter-rot',
                legend_height = unit(3, "cm"),
                adjust = c("right", "top")
              ))

draw(ht,annotation_legend_side = "right",row_dend_side = "left")

