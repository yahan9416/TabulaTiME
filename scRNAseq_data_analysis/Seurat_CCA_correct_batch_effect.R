#' Correct batch effect by CCA and return a seurat object.
#' @param H5_files Path of expression matrix of h5 file
#' @param Seurat_obj If you have built a Seurat object, and you can use it as input
#' @param Batch_info A vector of batch information corresponding to cell barcode in expmat or seurat object
#' @param nfeatures How many feature will be used for find anchors and appeared in the integrated object.
#' @param Gene_list A vector of a rds file including genes which you want to keep in the downstream analysis
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)

Run_CCA_Specific_gene<-function(H5_files=NULL,Seurat_obj=NULL,Batch_info,nfeatures=5000,Gene_list=NULL,dims.use=1:30){
  runpca.agrs = list(npcs = 100)
  findneighbors.args = list() 
  findclusters.args = list()
  
  if(length(H5_files) >0){
    # read data
    expr = Read10X_h5(H5_files)
    seuratObj <- CreateSeuratObject(counts = expr, project = "Seurat_obj", min.cells = 10, min.features = 500)
    seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
    seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj$RNA))
    seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj),npcs = 100)
    
    Seurat_obj=seuratObj
  }
  Seurat_obj@meta.data$Batch=Batch_info
  data.list <- SplitObject(RNA, split.by = "batch")
  if(length(Gene_list) == 1){
    Gene_list=readRDS(Gene_list)
  }
  
  for(i in 1:length(data.list)){
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = (nfeatures - length(Gene_list)), verbose = FALSE)
    VariableFeatures(data.list[[i]])<-c(VariableFeatures(data.list[[i]]),interaction_protein)
  }
  
  anchors <- FindIntegrationAnchors(object.list = data.list,dims = dims.use, anchor.features = nfeatures)
  RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
  RNA.integrated@project.name <- paste0(RNA@project.name,'_CCA')
  DefaultAssay(RNA.integrated) <- "integrated"
  
  RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
  RNA.integrated <- fastDoCall("RunPCA", c(object = RNA.integrated, runpca.agrs))
  
  p = ElbowPlot(object = RNA.integrated, ndims = RNA.integrated@commands$RunPCA.integrated@params$npcs)
  ggsave(file.path(paste0(RNA.integrated@project.name, "_PCElbowPlot.png")), p,  width=10, height=4)
  
  RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", dims = dims.use, ...)
  RNA.integrated <- fastDoCall("FindNeighbors", c(object = RNA.integrated, reduction = "pca", dims = dims.use, findneighbors.args))
  RNA.integrated <- fastDoCall("FindClusters", c(object = RNA.integrated, resolution = cluster.res, findclusters.args))
  return(RNA.integrated)
}
