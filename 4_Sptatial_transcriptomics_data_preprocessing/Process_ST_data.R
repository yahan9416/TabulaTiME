#' Processing ST dataset
#' @param ST_Feature_bc_matrix Path of 10X generated ST datasets
#' @param num_cores Number of core used
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)

args <- commandArgs(trailingOnly = TRUE)
ST_Feature_bc_matrix<- args[1]
num_cores<- args[2]

ST_seurat=Load10X_Spatial(ST_Feature_bc_matrix)
ST_seurat <- PercentageFeatureSet(ST_seurat, "^MT-", col.name = "percent_mito")
ST_seurat = ST_seurat[, ST_seurat$nFeature_Spatial > 200 & ST_seurat$percent_mito < 25 ]
ST_seurat <- ST_seurat[!grepl("^mt-", rownames(ST_seurat)), ]
####quality control
pdf("ST_Seurat_QC.pdf",width = 12,height = 12)
SpatialFeaturePlot(ST_seurat, features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"))
dev.off()


######Marker gene Visualization
DefaultAssay(ST_seurat)<-"Spatial"
SpatialFeaturePlot(object = ST_seurat, features = c("COL1A1","CTHRC1","SFRP1","KRT19","PTPRC","CCL5","SAA1"), ncol = num_cores)
