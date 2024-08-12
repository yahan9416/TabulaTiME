#Fibroblast as example

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan(multisession, workers = 12)
options(future.globals.maxSize = 10*1024^4)

AllMetacell_SeuratObj=readRDS("./PanCancer_Minicluster30_Seurat_obj.rds")
Fibroblast_meta=readRDS("./SeuratObj_meta_info_addTeatmentResp.rds")
Fibroblast_cell=intersect(rownames(SeuratObj@meta.data),rownames(Fibroblast_meta)[which(Fibroblast_meta$curated_annota %in% c("Fibroblast","Myofibroblast"))])
SeuratObj=subset(SeuratObj,cells=Fibroblast_cell)
Dataset_info=read.table("./TISCH_V1V2_10X_data_patient_cell1.txt",header = TRUE,sep="\t")
SeuratObj@meta.data$Tissue=Dataset_info$Tissue[match(SeuratObj@meta.data$Dataset,Dataset_info[,1])]

RNA=SeuratObj
batch=SeuratObj@meta.data$Tissue
nfeatures = 5000
dims.use = 1:30
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
Metabolic_gene=readRDS("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Lineage_seperate_analysis/KEGG_enzyme_metabolic_gene_list.rds")
Metabolic_gene=unique(Metabolic_gene)
data.list <- SplitObject(RNA, split.by = "batch")
for(i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 3333, verbose = FALSE)
  VariableFeatures(data.list[[i]])<-c(VariableFeatures(data.list[[i]]),Metabolic_gene)
}
anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  dims = dims.use, anchor.features = nfeatures)

RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
RNA.integrated@project.name <- "Minicluster30_Fibroblast_CCA_addMeta"
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
