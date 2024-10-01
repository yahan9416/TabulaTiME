#' Utilize copykat to predict the presence of malignant cells
#' @param Seurat_obj_path Path of seurat object
#' @param num_cores Number of core used
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj_path <- args[1]
num_cores <-args[2]


library(copykat)
library(Seurat)

Epi_Seurat=readRDS(Seurat_obj_path)
#DefaultAssay(Epi_Seurat)<-"integrated"
expmat=as.matrix(Epi_Seurat@assays$RNA@counts)
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=num_cores)

predict=copykat.test$prediction[,2]