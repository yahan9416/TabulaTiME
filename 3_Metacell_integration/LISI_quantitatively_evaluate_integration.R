library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(lisi)
Seurat_obj
umap_df = Seurat_obj$RNA@reductions$umap@cell.embeddings
#umap_celltype=Seurat_obj$RNA@meta.data$Dataset
res <- compute_lisi(umap_df, Seurat_obj$RNA@meta.data, c('Dataset'))
