library(cluster)
library(Seurat)
library(MAESTRO)
library(clustree)
SeuratObj=readRDS("./Merge_minicluster30_Fibroblast_Seurat_CCA_addMeta.rds")

cluster_label=SeuratObj@meta.data$seurat_clusters
cluster_label=as.numeric(cluster_label)
rd=SeuratObj@reductions$umap@cell.embeddings
si_umap_can <- silhouette(cluster_label, dist(rd, "canberra"))
si_umap_eu <- silhouette(cluster_label, dist(rd, "euclidean"))
si_umap_max <- silhouette(cluster_label, dist(rd, "maximum"))

rd=SeuratObj@reductions$pca@cell.embeddings[,1:10]
si_pca_can <- silhouette(cluster_label, dist(rd, "canberra"))
si_pca_eu <- silhouette(cluster_label, dist(rd, "euclidean"))
si_pca_max <- silhouette(cluster_label, dist(rd, "maximum"))