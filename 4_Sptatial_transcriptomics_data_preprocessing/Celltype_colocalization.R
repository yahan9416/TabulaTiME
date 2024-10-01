#' Evaluating the cell type colocalozation in ST dataset 
#' The positive correlation means the cell types trend to appear in one spot
#' @param Celltype1_DEG_gene_list_path Path of DEG files
#' @param Celltype2_DEG_gene_list_path Path of DEG files
#' @param Deconvolution_result_path Path of the deconvolution result rom STRIDE
#' @param Seurat_obj_path Path of Seurat object
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Celltype1_DEG_gene_list_path<-args[1]
Celltype2_DEG_gene_list_path<-args[2]
Deconvolution_result_path<-args[3]
Seurat_obj_path<-args[4]

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 10000 * 1024^4)

DEG_gene_list=readRDS(Celltype1_DEG_gene_list_path)
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


DEG_gene_list=readRDS(Celltype2_DEG_gene_list_path)
DEG_gene_list=droplevels(DEG_gene_list)
DEG_gene_list=DEG_gene_list[ order(DEG_gene_list$avg_log2FC,decreasing = TRUE),]
Mye_DEG_gene_list<<-list()
process_DEG_list_form<-function(temp_cluster){
  temp_DEG=DEG_gene_list$gene[which(DEG_gene_list$cluster == temp_cluster)]
  temp_DEG=setdiff(temp_DEG,unlist(Fibro_DEG_gene_list))
  if(length(temp_DEG) > 50){
    Mye_DEG_gene_list<<-c(Mye_DEG_gene_list,list(temp_DEG[1:50]))
  }else{
    Mye_DEG_gene_list<<-c(Mye_DEG_gene_list,list(temp_DEG))}
}
result=apply(matrix(unique(as.vector(unlist(DEG_gene_list$cluster)))),1,process_DEG_list_form)
names(Mye_DEG_gene_list)=unique(as.vector(unlist(DEG_gene_list$cluster)))


PLC_Myeloid_CTHRC1_cor<<-NULL
  HCCIL_fra=readRDS(Deconvolution_result_path) #load the deconvolution result rom STRIDE
  Seurat_obj=readRDS(Seurat_obj_path) #load the seurat object
  
  
  co_existFM=intersect(which(HCCIL_fra$Fibroblasts > 0.15),which(HCCIL_fra$Mono.Macro > 0.15))
  if(length(co_existFM) < 10){
    co_existFM=intersect(which(HCCIL_fra$Fibroblasts > 0.1),which(HCCIL_fra$Mono.Macro > 0.1))
  }
  
  Seurat_obj=AddModuleScore(
    object = Seurat_obj,
    features = c(Fibro_DEG_gene_list,list(Mye_DEG_gene_list$Macro_SLPI)),
    nbin=20,
    assay	="Spatial",
    ctrl = 5,
    name = c(names(Fibro_DEG_gene_list),"Macro_SLPI"))
  
  colnames(Seurat_obj@meta.data)[5:16]=c(names(Fibro_DEG_gene_list),"Macro_SLPI")
  batch_calculate_correlation<-function(temp_index){
    print(temp_index)
    cor_result=unlist(cor.test(Seurat_obj@meta.data[,"Fibro_CTHRC1"],Seurat_obj@meta.data[,temp_index]))
       temp_cor=c(temp_file,"Macro_SLPI",colnames(Seurat_obj@meta.data)[temp_index],cor_result[4],cor_result[3])
       PLC_Myeloid_CTHRC1_cor<<-rbind(PLC_Myeloid_CTHRC1_cor,temp_cor)
    #return(c("Fibro_CTHRC",colnames(Seurat_obj@meta.data)[temp_index],cor_result[4],cor_result[3]))
  }
  result=apply(matrix(5:15),1,batch_calculate_correlation)  
  

colnames(PLC_Myeloid_CTHRC1_cor)=c("Sample","Fibro","Myeloid","Correlation","Pvalue")

#The positive correlation means the cell types trend to appear in one spot