#' Calculating the LISI score to evaluate the performance of batch effect correction.
#' @param CCA_seurat_obj Path of CCA corrected Seurat obj
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(lisi)

args <- commandArgs(trailingOnly = TRUE)
CCA_seurat_obj<-args[1]

Seurat_obj=readRDS(CCA_seurat_obj)
umap_df = Seurat_obj$RNA@reductions$umap@cell.embeddings
#umap_celltype=Seurat_obj$RNA@meta.data$Dataset
res <- compute_lisi(umap_df, Seurat_obj$RNA@meta.data, c('Dataset'))
