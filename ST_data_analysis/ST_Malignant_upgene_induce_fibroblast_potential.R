#which malignant up-regulated genes may induce the activation of fibroblasts

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 10000 * 1024^4)

#After generate the Seurat object for each sample
DEG_gene_list=readRDS("./Fibroblast/Basic_analysis/Fibroblast_celltype_DEG.rds")
DEG_gene_list=droplevels(DEG_gene_list)
DEG_gene_list=DEG_gene_list[ order(DEG_gene_list$avg_log2FC,decreasing = TRUE),]
Fibro_DEG_gene_list<<-list()
process_DEG_list_form<-function(temp_cluster){
  temp_DEG=DEG_gene_list$gene[which(DEG_gene_list$cluster == temp_cluster)]
  if(length(temp_DEG) > 50){
    Fibro_DEG_gene_list<<-c(Fibro_DEG_gene_list,list(temp_DEG[1:50]))
  }else{
    Fibro_DEG_gene_list<<-c(Fibro_DEG_gene_list,list(temp_DEG))}
}
result=apply(matrix(unique(DEG_gene_list$cluster)),1,process_DEG_list_form)
names(Fibro_DEG_gene_list)=unique(DEG_gene_list$cluster)

Seurat_obj #load the Seurat object of your want to estimated
CopyKat=readRDS("CopyKat.rds"))
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
