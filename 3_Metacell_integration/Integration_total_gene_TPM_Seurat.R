#' Integrated gene expression matrix
#' @param Metacell_Expression_obj Path of integrated metacell expression matrix
#' @param Names_Seuratobj Name of integrated Seurat object
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future) 
future::plan(multisession, workers = 1)
options(future.globals.maxSize = 10*1024^4)


args <- commandArgs(trailingOnly = TRUE)
Metacell_Expression_obj<-args[1]
Names_Seuratobj<-args[2]

Minicluster_TPM=readRDS(Metacell_Expression_obj)
SeuratObj <- CreateSeuratObject(Minicluster_TPM, project = Names_Seuratobj, 
                                min.cells = 10, min.features = 500)
mito.genes <- grep("^MT-", rownames(GetAssayData(object = SeuratObj)), 
                   value = TRUE)
ercc.genes <- grep("^ERCC", rownames(GetAssayData(object = SeuratObj)), 
                   value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = SeuratObj)[mito.genes, 
                                                                 ])/Matrix::colSums(GetAssayData(object = SeuratObj))
percent.ercc <- Matrix::colSums(GetAssayData(object = SeuratObj)[ercc.genes, 
                                                                 ])/Matrix::colSums(GetAssayData(object = SeuratObj))
SeuratObj$percent.mito <- percent.mito
SeuratObj$percent.ercc <- percent.ercc
p1 = VlnPlot(SeuratObj, c("percent.mito", "percent.ercc"), 
             ncol = 2)
ggsave(file.path( paste0(SeuratObj@project.name, 
                         ".spikein.png")), p1, width = 6, height = 4.5)
SeuratObj <- subset(x = SeuratObj, subset = percent.mito < 
                      0.15)
SeuratObj <- subset(x = SeuratObj, subset = percent.ercc < 
                      0.05)
vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ercc")
message("Normalization and identify variable genes ...")

SeuratObj <- NormalizeData(object = SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

#RowsNA_gene<-names(which(rowSums(is.na(SeuratObj@assays@counts))>0))  

SeuratObj <- FindVariableFeatures(object = SeuratObj, selection.method = "vst",   nfeatures = 3000)
SeuratObj <- ScaleData(object = SeuratObj)
message("PCA analysis ...")
library(Gmisc)
SeuratObj <- fastDoCall("RunPCA", c(object = SeuratObj, features = VariableFeatures(SeuratObj), 
                                    runpca.agrs=list()))
#p2 = ElbowPlot(object = SeuratObj, ndims = SeuratObj@commands$RunPCA.RNA@params$npcs)
#ggsave(file.path(outdir, paste0(SeuratObj@project.name, "_PCElbowPlot.png")), 
#    p2, width = 10, height = 4)
message("UMAP analysis ...")
SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", 
                     dims = 1:50)
SeuratObj <- fastDoCall("FindNeighbors", c(object = SeuratObj, 
                                           reduction = "pca", dims = 1:50, findneighbors.args=list()))


saveRDS(SeuratObj,file="PanCancer_Metacell_Seurat_obj.rds")  