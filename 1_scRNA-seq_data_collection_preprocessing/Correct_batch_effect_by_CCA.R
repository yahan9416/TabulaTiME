#' Correct the batch effect by CCA and return a seurat object.
#' @param Seurat_obj_path The path of Seurat object
#' @param Batch_info A column name of meta.data, including the batch information
#' @param Num_features Number of features used for integration 
#' @param Num_features Number of dims used for integration
#' @param ncores Number of cores used 
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj_path <- args[1]
Num_features <- args[2]
Num_dimused <- args[3]
Batch_info <- args[4] #the column of meta.data represent batch
ncores <- args[5]

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 10*1024^4)


seuratObj=readRDS(Seurat_obj_path)
batch=seuratObj@meta.data[,match(Batch_info,colnames(seuratObj@meta.data))]
nfeatures = Num_features
dims.use = 1:Num_dimused
cluster.res =1
only.pos = FALSE
genes.test.use = "presto"
genes.cutoff = 1e-05
genes.pct = 0.1
genes.logfc = 0.25
runpca.agrs = list()
findneighbors.args = list()
findclusters.args = list()
runpca.agrs = list(npcs = 50)
only.pos = TRUE
seuratObj@meta.data$batch <- batch
data.list <- SplitObject(seuratObj, split.by = "batch")
for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], 
                                         selection.method = "vst", nfeatures = nfeatures, 
                                         verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  dims = dims.use, anchor.features = nfeatures)
seuratObj.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
seuratObj.integrated@project.name <- paste0(seuratObj@project.name, "_CCA")
DefaultAssay(seuratObj.integrated) <- "integrated"
seuratObj.integrated <- ScaleData(seuratObj.integrated, verbose = FALSE)
seuratObj.integrated <- fastDoCall("RunPCA", c(object = seuratObj.integrated, 
                                         runpca.agrs))
p = ElbowPlot(object = seuratObj.integrated, ndims = seuratObj.integrated@commands$RunPCA.integrated@params$npcs)
ggsave(file.path(paste0(seuratObj.integrated@project.name, "_PCElbowPlot.png")), 
       p, width = 10, height = 4)
seuratObj.integrated <- RunUMAP(object = seuratObj.integrated, reduction = "pca", 
                          dims = dims.use)
seuratObj.integrated <- fastDoCall("FindNeighbors", c(object = seuratObj.integrated, 
                                                reduction = "pca", dims = dims.use, findneighbors.args))
seuratObj.integrated <- fastDoCall("FindClusters", c(object = seuratObj.integrated, 
                                               resolution = cluster.res, findclusters.args))
seuratObj.integrated@project.name="GSE_Seurat_CCA"

p = DimPlot(object = seuratObj.integrated, label = TRUE, pt.size = 0.2)

ggsave(file.path(paste0(seuratObj.integrated@project.name, "_cluster.png")), 
       p, width = 5, height = 4)
p = DimPlot(object = seuratObj.integrated, group = "batch", label = TRUE, 
            pt.size = 0.2)
ggsave(file.path(paste0(seuratObj.integrated@project.name, "_batch.png")), 
       p, width = 5.5, height = 4)
cluster.genes <- FindAllMarkersMAESTRO(object = seuratObj.integrated, 
                                       only.pos = only.pos, min.pct = genes.pct, test.use = genes.test.use, 
                                       logfc.threshold = genes.logfc)
cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 
                                 genes.cutoff, ]
write.table(cluster.genes, paste0(seuratObj.integrated@project.name, 
                                  "_DiffGenes.tsv"), quote = F, sep = "\t")
saveRDS(seuratObj.integrated,file="Seuratobjt_CCA.rds")