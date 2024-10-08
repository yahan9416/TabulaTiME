#' Clustree were used to select the optimal cluster resolution
#' @param Seurat_obj Path of Seurat object
#' @author Ya Han

library(Seurat)
library(MAESTRO)
library(clustree)
args <- commandArgs(trailingOnly = TRUE)
Seurat_obj <- args[1]


SeuratObj=readRDS(Seurat_obj)
SeuratObj <- FindClusters(
  object = SeuratObj,
  resolution = c(seq(0.1,1.5,0.1))
)

pdf("Clustree_resolution_marker_Gene_exoression.pdf",height = 15,width = 15)
clustree(SeuratObj, prefix = "integrated_snn_res.",node_colour = "Marker_gene",
         node_colour_aggr = "mean")
dev.off()
