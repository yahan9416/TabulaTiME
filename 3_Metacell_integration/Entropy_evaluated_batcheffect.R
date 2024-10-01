#' Calculating the Entropy score to evaluate the performance of batch effect correction.
#' @param Seurat_obj Path of integrated seurat object
#' @author Ya Han

library(mclust)
library(rjson)
library(dbscan)
library(MAESTRO)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj <- args[1]

Sc_V2_Seurat=readRDS(Seurat_obj)

set.seed(0)
k = 30
knn <- kNN(Sc_V2_Seurat@reductions[["umap"]]@cell.embeddings,k = k)
neighbor_matrix <- knn$id
#fetch a vector of batch conditions: ["patient1","patient1","patient2"...]
Vbatch <- Sc_V2_Seurat@meta.data$Patient
#batch compositions: ["patient1","patient2","patient3"...]
batch <- unique(Vbatch)
Entropies_Sc_V2_Seurat = sapply(1:ncol(Sc_V2_Seurat), function(i){
  sum=0
  cell_neighbors <- Vbatch[neighbor_matrix[i,]]
  batch_pct = sapply(1:length(batch), function(j){
    Nbatch <- length(which(cell_neighbors==batch[j]))
    if (Nbatch == 0){
      return(0)
    } else {
      return((Nbatch/k)*log2(Nbatch/k))
    }
  })
  return(-sum(batch_pct))
})

saveRDS(Entropies_Sc_V2_Seurat,file="Entropies_Sc_V2_Seurat.rds")