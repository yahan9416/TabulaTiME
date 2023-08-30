#' DoubletFinder is tools to identify the doublets for each single cell RNAseq data
#' @param H5_files Path of h5 file folder 
#' @author Ya Han

library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(future)

Run_DoubletFinder<-function(H5_files){
  files=list.files(H5_files,recursive=T)
  files=files[grep(".h5",files)]
  Batch_identify_Doublet_for_each_sample<-function(x){
    expr = Read10X_h5(paste0("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Data/",x))
    Temp_Seurat <- CreateSeuratObject(expr)
    Temp_Seurat <- NormalizeData(Temp_Seurat)
    Temp_Seurat <- FindVariableFeatures(Temp_Seurat, selection.method = "vst", nfeatures = 2000)
    Temp_Seurat <- ScaleData(Temp_Seurat)
    Temp_Seurat <- RunPCA(Temp_Seurat)
    Temp_Seurat <- RunUMAP(Temp_Seurat, dims = 1:10)
    ## pK Identification (no ground-truth)
    sweep.res.list_kidney <- paramSweep_v3(Temp_Seurat, PCs = 1:10, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    Temp_Seurat <- FindNeighbors(Temp_Seurat, dims = 1:10)
    Temp_Seurat <- FindClusters(Temp_Seurat, resolution = 0.5)
    annotations <- Temp_Seurat@meta.data$seurat_clusters 
    homotypic.prop <- modelHomotypic(annotations)           ## ex: 
    nExp_poi <- round(0.075*dim(Temp_Seurat@meta.data)[1])  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    Temp_Seurat <- doubletFinder_v3(Temp_Seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    return(cbind(rownames(Temp_Seurat@meta.data),Temp_Seurat@meta.data[,7]))
  }
  result=apply(matrix(files),1,Batch_identify_Doublet_for_each_sample)
  sample=lapply(strsplit(files,"/"),function(x) x[1])
  names(result)=unlist(sample)
  return(result)
}