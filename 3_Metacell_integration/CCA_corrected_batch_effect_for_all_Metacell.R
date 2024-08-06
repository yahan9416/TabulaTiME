library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 10*1024^4)
seuratObj=readRDS("。/PanCancer_Metacell_30_Seurat_obj.rds")

seuratObj@meta.data$Dataset=unlist(lapply(strsplit(rownames(seuratObj@meta.data),"@"),function(x) x[1]))
seuratObj@meta.data$Cancer_type=unlist(lapply(strsplit(rownames(seuratObj@meta.data),"_"),function(x) x[1]))
Dataset_info=read.table("TISCH_V1V2_10X_data_patient_cell.txt",header = TRUE,sep="\t")
seuratObj@meta.data$Tissue=Dataset_info$Tissue[match(seuratObj@meta.data$Dataset,Dataset_info[,1])]


RNA=seuratObj
batch=seuratObj@meta.data$Tissue
nfeatures = 3000
dims.use = 1:50
cluster.res =1
only.pos = FALSE
genes.test.use = "presto"
genes.cutoff = 1e-05
genes.pct = 0.1
genes.logfc = 0.25
runpca.agrs = list()
findneighbors.args = list()
findclusters.args = list()
runpca.agrs = list(npcs = 75)
only.pos = TRUE
RNA@meta.data$batch <- batch
data.list <- SplitObject(RNA, split.by = "batch")
for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], 
                                         selection.method = "vst", nfeatures = nfeatures, 
                                         verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  dims = dims.use, anchor.features = nfeatures)
#saveRDS(anchors,file="Merge_minicluster30_Seurat_CCA_anchors.rds")
#anchors=readRDS("Merge_minicluster30_Seurat_CCA_anchors.rds")
RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
RNA.integrated@project.name <- paste0(RNA@project.name, "_CCA")
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
RNA.integrated@project.name="Merge_minicluster30_Seurat_CCA"

p = DimPlot(object = RNA.integrated, label = TRUE, pt.size = 0.2)

ggsave(file.path(paste0(RNA.integrated@project.name, "_cluster.png")), 
       p, width = 5, height = 4)
p = DimPlot(object = RNA.integrated, group = "batch", label = TRUE, 
            pt.size = 0.2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_batch.png")), 
       p, width = 5.5, height = 4)
cluster.genes <- FindAllMarkersMAESTRO(object = RNA.integrated, 
                                       only.pos = only.pos, min.pct = genes.pct, test.use = genes.test.use, 
                                       logfc.threshold = genes.logfc)
cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 
                                 genes.cutoff, ]
write.table(cluster.genes, paste0(RNA.integrated@project.name, 
                                  "_DiffGenes.tsv"), quote = F, sep = "\t")

DefaultAssay(RNA.integrated)<-"RNA"
VisualizeUmap(SeuratObj = RNA.integrated, type = "RNA", genes = c("PTPRC", "CD8A", "CD8B", "CD3D","CD4","CD14","CD16","CD1C","LAMP3","CD79A","JACHIN","VWF","COL1A1","DCN","EPCAM","MYLK","PECAM1","MKI67","TOP2A"), ncol = 6,
              width = 30, height = 12, name = paste0(RNA.integrated@project.name, "_immunity"))
saveRDS(RNA.integrated,file="Merge_minicluster30_Seurat_CCA.rds")