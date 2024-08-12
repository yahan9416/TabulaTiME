library(cluster)
library(Seurat)
library(MAESTRO)
library(clustree)
args <- commandArgs(trailingOnly = TRUE)
CCA_seurat_obj<-args[1]

SeuratObj=readRDS(CCA_seurat_obj)

cluster_label=SeuratObj@meta.data$seurat_clusters
cluster_label=as.numeric(cluster_label)

rd=SeuratObj@reductions$pca@cell.embeddings[,1:30]
si_pca_can <- silhouette(cluster_label, dist(rd, "canberra"))
si_pca_eu <- silhouette(cluster_label, dist(rd, "euclidean"))
si_pca_max <- silhouette(cluster_label, dist(rd, "maximum"))