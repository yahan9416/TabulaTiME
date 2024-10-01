#' Correct batch effect by Harmony and return a seurat object.
#' @param H5_files Path of expression matrix of h5 file
#' @param Seurat_obj If you have built a Seurat object, and you can use it as input
#' @param Batch_info A vector of batch information corresponding to cell barcode in expmat or seurat object
#' @author Ya Han

library(harmony)
library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)

args <- commandArgs(trailingOnly = TRUE)
H5_files<-args[1]
Seurat_obj<-args[2]
Batch_info<-args[3]

Run_Harmony<-function(Seurat_obj= NULL,Batch_info){
  if(length(H5_files) >0){
    # read data
    seuratObj@meta.data$Batch=Batch_info
    
    seuratObj <- RunHarmony(seuratObj, "Batch")
    seuratObj=RunUMAP(seuratObj,reduction = "harmony", dims = 1:75)
    seuratObj <- FindNeighbors(seuratObj, dims = 1:75)
    seuratObj <- FindClusters(seuratObj, resolution = 1)
    
    return(seuratObj)
  }
  
  if(length(Seurat_obj) >0){
    Seurat_obj@meta.data$Batch=Batch_info
    Seurat_obj <- RunHarmony(Seurat_obj, "Batch")
    Seurat_obj=RunUMAP(Seurat_obj,reduction = "harmony", dims = 1:75)
    Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:75)
    Seurat_obj <- FindClusters(Seurat_obj, resolution = 1)
    return(Seurat_obj)
  }
  
}
