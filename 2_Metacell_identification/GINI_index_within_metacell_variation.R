#' Calculating the GINI index to evaluate the optimal number of cells within each metacell.
#' @param Metacell_info Path of metacell information
#' @param Seurat_obj_path Path of integrated seurat object
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(lawstat)


args <- commandArgs(trailingOnly = TRUE)
Metacell_info <- args[1]
Seurat_obj_path<-args[2]

load(Metacell_info)
Seurat_obj=readRDS(Seurat_obj_path)
Variation_withincluster<<-NULL

  #Now we need to generate dataframe
  cell_minicluster<<-NULL
  result=apply(matrix(1:length(colname_source_cluster)),1,function(y){ 
    cell_minicluster<<-rbind(cell_minicluster,cbind(Cell_minicluster_list[[y]],colname_source_cluster[y]))
  })
  Seurat_obj@meta.data$minicluster_id="0"
  Seurat_obj@meta.data$minicluster_id[match(cell_minicluster[,1],rownames(Seurat_obj@meta.data))]=cell_minicluster[,2]
  DefaultAssay(Seurat_obj$RNA)<-"RNA"
  count=GetAssayData(Seurat_obj$RNA,slot ="counts")
  cal_TPM <- RNACountToTPM(count, idType = "SYMBOL", organism = "GRCh38")
  log_TPM=log2((cal_TPM/10)+1)
  log_TPM=as.matrix(log_TPM)
  
  #calcluated Gini-index
  mini_cluster_gini<-function(temp_mini){
    temp_tpm=log_TPM[,which(Seurat_obj@meta.data$minicluster_id == temp_mini)]
    marker_gini <- mean(na.omit(apply(temp_tpm,1,function(x) gini.index(x)$statistic)))
    return(marker_gini)
  }
  result=apply(matrix(setdiff(unique(unlist(Seurat_obj@meta.data$minicluster_id)),"0")[1:3]),1,mini_cluster_gini)
  
saveRDS(result,"Metacell_corresponding_Gini_index.rds")
