library(msigdbr)

msigdbr_show_species()

### C2
c2 = msigdbr(species = "Rattus norvegicus", category = "C2")
c2
c2 <- c2 %>% group_by(gs_name) %>% 
  summarise(gene=paste0(gene_symbol, collapse = '\t'))

c2 <- paste0(c2$gs_name, '\t', '> ', c2$gs_name, '\t', c2$gene)
write.table(c2, file='data/msigdbr/rat/c2.all.rat.symbols.gmt', row.names=F, col.names = F, quote = F)


### C2 KEGG
c2.kegg = msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = 'CP.KEGG')
c2.kegg
c2 <- c2 %>% group_by(gs_name) %>% 
  summarise(gene=paste0(gene_symbol, collapse = '\t'))

c2 <- paste0(c2$gs_name, '\t', '> ', c2$gs_name, '\t', c2$gene)
write.table(c2, file='data/msigdbr/rat/c2.all.rat.symbols.gmt', row.names=F, col.names = F, quote = F)


### C5 BP
c5.bp = msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = 'BP')
c5.bp

c5.bp <- c5.bp %>% group_by(gs_name) %>% 
  summarise(gene=paste0(gene_symbol, collapse = '\t'))

c5.bp <- paste0(c5.bp$gs_name, '\t', '> ', c5.bp$gs_name, '\t', c5.bp$gene)
write.table(c5.bp, file='data/msigdbr/rat/c5.bp.rat.symbols.gmt', row.names=F, col.names = F, quote = F)


### HALLMARK
hallmark = msigdbr(species = "Rattus norvegicus", category = "H")
hallmark

hallmark <- hallmark %>% group_by(gs_name) %>% 
  summarise(gene=paste0(gene_symbol, collapse = '\t'))

hallmark <- paste0(hallmark$gs_name, '\t', '> ', hallmark$gs_name, '\t', hallmark$gene)
write.table(hallmark, file='data/msigdbr/rat/h.all.rat.symbols.gmt', row.names=F, col.names = F, quote = F)



