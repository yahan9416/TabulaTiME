#' Adjusted Rand Index (ARI) was used to evaluate the performance of batch effect correction.
#' @param CCA_seurat_obj Path of CCA corrected seurat object
#' @author Ya Han

library(mclust)
library(MAESTRO)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
CCA_seurat_obj<-args[1]

Seurat_obj=readRDS(CCA_seurat_obj)
set.seed(0)
#1-ARIbatch
cluster_label=Seurat_obj@meta.data$seurat_clusters
cluster_blebe=Seurat_obj@meta.data$Celltype
batch_random_20<-function(temp_index){
  temp_index=sample(1:length(cluster_batch),5000)
  ari_batch <- adjustedRandIndex(cluster_batch[temp_index], cluster_label[temp_index])
  return(ari_batch)
}
result=apply(matrix(1:20),1,batch_random_20)
saveRDS(result,"ARI_batch.rds")
#ARIcelltype
cluster_label=Seurat_obj@meta.data$seurat_clusters
cluster_batch=Seurat_obj@meta.data$annotation
batch_random_20<-function(temp_index){
  temp_index=sample(1:length(cluster_batch),5000)
  ari_celltype <- adjustedRandIndex(cluster_batch[temp_index], cluster_label[temp_index])
  return(ari_celltype)
}
result=apply(matrix(1:20),1,batch_random_20)
saveRDS(result,"ARI_celltypes.rds")