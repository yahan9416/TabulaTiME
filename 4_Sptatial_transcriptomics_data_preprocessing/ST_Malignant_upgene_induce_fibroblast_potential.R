#' Infer the correlation between up-regulated malignant genes and the activation of fibroblasts, aiming to identify the genes that may induce fibroblast activation.
#' @param DEG_gene_path Path of cell type specific DEG file
#' @param Seurat_obj_path Path of Seurat object
#' @param CopyKat_result_path Path of copykat predict result
#' @param num_cores Number of core used
#' @author Ya Han

##CopyKAT######
args <- commandArgs(trailingOnly = TRUE)
DEG_gene_path<- args[1]
Seurat_obj_path<- args[2]
CopyKat_result_path<- args[3]
num_cores<- args[4]



library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = num_cores)
options(future.globals.maxSize = 10000 * 1024^4)


#After generate the Seurat object for each sample
Fibro_DEG_gene_list=readRDS(DEG_gene_path)

Seurat_obj=readRDS(Seurat_obj_path) #load the Seurat object of your want to estimated
CopyKat=readRDS(CopyKat_result_path)
Seurat_obj=ScaleData(Seurat_obj,verbose = FALSE)
Seurat_obj=FindVariableFeatures(Seurat_obj,nfeatures=3000)
Seurat_obj <- RunPCA(Seurat_obj, verbose = FALSE)
Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:30)
Seurat_obj <- FindClusters(Seurat_obj, verbose = FALSE)
Seurat_obj <- RunUMAP(Seurat_obj, dims = 1:30)
de_markers <- FindAllMarkers(Seurat_obj)

Seurat_obj=AddModuleScore(
  object = Seurat_obj,
  features = Fibro_DEG_gene_list,
  nbin=20,
  assay	="Spatial",
  ctrl = 5,
  name = names(Fibro_DEG_gene_list))

colnames(Seurat_obj@meta.data)[(dim(Seurat_obj@meta.data)[2]-10):dim(Seurat_obj@meta.data)[2]]=names(Fibro_DEG_gene_list)

expmat=GetAssayData(Seurat_obj)
tumor_cell=rownames(CopyKat$prediction)[grep("aneuploid",CopyKat$prediction[,2])]
tumor_fra_copykat=length(tumor_cell)/dim(Seurat_obj)[2]
Seurat_obj@meta.data$tumor_cell="Normal"
Seurat_obj@meta.data$tumor_cell[which(rownames(Seurat_obj@meta.data) %in% tumor_cell)]="Tumor"
Tumor_cluster=table(Seurat_obj@meta.data$seurat_clusters,Seurat_obj@meta.data$tumor_cell)
Tumor_cluster=rownames(Tumor_cluster)[which(Tumor_cluster[,2] > Tumor_cluster[,1])]


Fibro_CTHRC1_cells=rownames(Seurat_obj@meta.data)[which(Seurat_obj@meta.data$Seurat_clusters  %in% c(4,11))]
loc_Tumor=GetTissueCoordinates(Seurat_obj, scale = NULL)
rownames(loc_Tumor)=rownames(Seurat_obj@meta.data)
tumor_cell_loc=loc_Tumor[which(rownames(loc_Tumor) %in% tumor_cell),]
Fibro_CTHRC1_loc=loc_Tumor[which(rownames(loc_Tumor) %in% tumor_cell),]


de_markers=de_markers[which(de_markers$avg_log2FC > 0),]
Mali_DEG=de_markers$gene[which(de_markers$cluster %in% Tumor_cluster)]

#Remove the cluster with higher expressed the marker gene of Fibro_CTHRC1: cluster 4 and 11
Mali_DEG=setdiff(Mali_DEG,de_markers$gene[which(de_markers$cluster %in%  c(4,11) )])
Mali_DEG=setdiff(Mali_DEG,Fibro_DEG_gene_list$FibroCTHRC1)

expmat=as.matrix(expmat)
expmat[which(is.na(expmat))]=0
gene_exp_mean_cthrc1=rowMeans(expmat[,match(rownames(Seurat_obj@meta.data)[which(Seurat_obj@meta.data$seurat_clusters %in% c(4,11))], colnames(expmat))])
CTHRC1_high=names(gene_exp_mean_cthrc1)[which(gene_exp_mean_cthrc1 > 1)]
Mali_DEG=setdiff(Mali_DEG,CTHRC1_high)

Malig_DEG_FibroCTHRC1_cor<<-NULL
Find_Mali_induce_correlation_Fibro_CTHRC1<-function(temp_gene){
  batch_calculate_distance_Fibro<-function(temp_tumor_cell){
    tumor_inedx=match(temp_tumor_cell,rownames(tumor_cell_loc))
    tumor_dis=mean(sqrt((tumor_cell_loc$imagerow[tumor_inedx]-Fibro_CTHRC1_loc$imagerow)^2+(tumor_cell_loc$imagecol[tumor_inedx]-Fibro_CTHRC1_loc$imagecol)^2))
  }
  Fibro_distance=apply(matrix(tumor_cell),1,batch_calculate_distance_Fibro)
  Mali_DEG_exp=expmat[match(temp_gene,rownames(expmat)),match(tumor_cell,colnames(expmat))]
  
  
  cor.test(Mali_DEG_exp,Fibro_distance)
  temp_result=c(temp_file,temp_gene,cor.test(Mali_DEG_exp,Fibro_distance))
  Malig_DEG_FibroCTHRC1_cor<<-rbind(Malig_DEG_FibroCTHRC1_cor,temp_result)
}
result=apply(matrix(unique(Mali_DEG)),1,Find_Mali_induce_correlation_Fibro_CTHRC1)
