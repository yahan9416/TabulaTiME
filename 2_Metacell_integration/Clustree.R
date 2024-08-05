library(Seurat)
library(MAESTRO)
library(clustree)
SeuratObj=readRDS("./Merge_minicluster30_Fibroblast_Seurat_CCA_addMeta.rds")
SeuratObj <- FindClusters(
  object = SeuratObj,
  resolution = c(seq(0.1,1.5,0.1))
)


Marker_gene="CTHRC1" #any intersted gene
pdf("Clustree_resolution_marker_Gene_exoression.pdf",height = 15,width = 15)
clustree(SeuratObj, prefix = "integrated_snn_res.",node_colour = "Marker_gene",
         node_colour_aggr = "mean")
dev.off()
