library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(lawstat)
Seurat_obj
minicluster_cellnumber#file list of minicluster_files
#minicluster_cellnumber=as.matrix(list.dirs("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Test_miniCluster_Performance/Minicluster"))
#for each number of cells in each mini-cluster
Variation_withincluster<<-NULL

batch_get_minicluster_distribution<-function(y){
  dataset_name=list.files(y,full.names = T)
  load(dataset_name[grep("CHOL_GSE125449",dataset_name)])
  
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
  rm(count)
  rm(cal_TPM)
  gc()
  #针对基因计算Gini-index
  mini_cluster_gini<-function(temp_mini){
    #match(rownames(Seurat_obj@meta.data)[which(Seurat_obj@meta.data$minicluster_id == temp_mini)],colnames(log_TPM))
    temp_tpm=log_TPM[,which(Seurat_obj@meta.data$minicluster_id == temp_mini)]
    marker_gini <- mean(na.omit(apply(temp_tpm,1,function(x) gini.index(x)$statistic)))
    return(marker_gini)
  }
  result=apply(matrix(setdiff(unique(unlist(Seurat_obj@meta.data$minicluster_id)),"0")[1:3]),1,mini_cluster_gini)
  setwd("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Test_miniCluster_Performance/Performance/Variation")
  saveRDS(result,file=paste0(Seurat_obj@project.name,basename(y),"_Gini_index.rds"))
  
  result=as.list(result)
  Variation_withincluster<<-c(Variation_withincluster,result)
}
result=apply(minicluster_cellnumber, 1, batch_get_minicluster_distribution)

