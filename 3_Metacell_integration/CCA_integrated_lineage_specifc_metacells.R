#' Correct batch effect by CCA and return a seurat object.
#' @param Seurat_obj If you have built a Seurat object, and you can use it as input
#' @param Batch_info A column name of meta.data, including the batch information
#' @param Nfeature Number of features used for integration 
#' @param N_dimused Number of dims used for integration
#' @param ncores Number of cores used 
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj<-args[1]
Batch_info<-args[2]
Nfeature<-args[3]
N_dimused<-args[4]
ncores<-args[5]

plan(multisession, workers = ncores)
options(future.globals.maxSize = 10*1024^4)

AllMetacell_SeuratObj=readRDS(Seurat_obj)

RNA=SeuratObj
batch=SeuratObj@meta.data[match(Batch_info,colnames(SeuratObj@meta.data))]
nfeatures = Nfeature
dims.use = 1:N_dimused
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
RNA@meta.data$batch <- batch
data.list <- SplitObject(RNA, split.by = "batch")
data.list <- SplitObject(RNA, split.by = "batch")
for(i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = Nfeature, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  dims = dims.use, anchor.features = nfeatures)

RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
RNA.integrated@project.name <- "Celllineage_CCA_addMeta"
DefaultAssay(RNA.integrated) <- "integrated"


RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
RNA.integrated <- fastDoCall("RunPCA", c(object = RNA.integrated, 
                                         runpca.agrs))
p = ElbowPlot(object = RNA.integrated, ndims = RNA.integrated@commands$RunPCA.integrated@params$npcs)
ggsave(file.path(paste0(RNA.integrated@project.name, "_PCElbowPlot.png")), 
       p, width = 10, height = 4)
RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", 
                          dims = dims.use)
RNA.integrated <- fastDoCall("FindNeighbors", c(object = RNA.integrated, 
                                                reduction = "pca", dims = dims.use, findneighbors.args))
RNA.integrated <- fastDoCall("FindClusters", c(object = RNA.integrated, 
                                               resolution = cluster.res, findclusters.args))
#RNA.integrated@project.name="Merge_minicluster30_Seurat_CCA"

p = DimPlot(object = RNA.integrated, label = TRUE, pt.size = 0.2)

umap_df = SeuratObj@reductions$pca
ggsave(file.path(paste0(RNA.integrated@project.name, "_cluster_new.png")), 
       p, width = 7, height = 4)
p = DimPlot(object = RNA.integrated, group = "batch", label = TRUE, 
            pt.size = 0.2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_batch.png")), 
       p, width = 8, height = 4)
p = DimPlot(object = RNA.integrated, group = "Source", label = TRUE, 
            pt.size = 0.2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_Source.png")), 
       p, width = 6, height = 4)

cluster.genes <- FindAllMarkersMAESTRO(object = RNA.integrated, 
                                       only.pos = only.pos, min.pct = genes.pct, test.use = genes.test.use, 
                                       logfc.threshold = genes.logfc)
cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 
                                 genes.cutoff, ]
write.table(cluster.genes, paste0(RNA.integrated@project.name, 
                                  "_DiffGenes_addMeta.tsv"), quote = F, sep = "\t")
