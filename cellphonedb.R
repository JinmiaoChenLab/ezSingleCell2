cellphone_for_seurat <- function(Seurat_obj){
counts <- as.data.frame(as.matrix(Seurat_obj@assays$RNA@data))

metadata <- data.frame(Cell = rownames(Seurat_obj@meta.data),cell_type = Seurat_obj@active.ident)
#metadata$Cell <- paste('d-pos_', metadata$Cell, sep = '')
write.table(counts,
file = 'counts.txt',
quote = F,
col.names = T,
row.names = T,
sep = '\t')
write.table(metadata,
file = 'metadata.txt',
quote = F,
col.names = T,
row.names = F,
sep = '\t')
system('cellphonedb method statistical_analysis --counts-data gene_name metadata.txt counts.txt --threads=50')
system('cellphonedb plot dot_plot')
system('cellphonedb plot heatmap_plot metadata.txt')
}