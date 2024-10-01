#' The Local Inverse Simpson Index (LISI) score was calculated to evaluate the performance of batch effect correction.
#' @param Integrated_seurat_obj Path of integrated seurat object
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(lisi)
args <- commandArgs(trailingOnly = TRUE)
Integrated_seurat_obj<-args[1]

Seurat_obj=readRDS(Integrated_seurat_obj)
umap_df = Seurat_obj$RNA@reductions$umap@cell.embeddings
#umap_celltype=Seurat_obj$RNA@meta.data$Dataset
res <- compute_lisi(umap_df, Seurat_obj$RNA@meta.data, c('Dataset'))
