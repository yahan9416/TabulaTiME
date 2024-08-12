library(Seurat)
library(MAESTRO)

args <- commandArgs(trailingOnly = TRUE)
Scrublet_result_path<- args[1]
Seurat_obj_path<- args[2]

#First collected all
setwd(Scrublet_result_path)
scrublet_result=list.files(Scrublet_result_path)
Scrublet_result<<-NULL
batch_identify_cellname<-function(temp_sample){
  sample_scrublet=read.table(temp_sample)
  sample_name=gsub("dbl_","",temp_sample)
  sample_name=gsub("matrix.mtx.txt","",sample_name)
  Dataset_sample=gsub("_count_","",sample_name)
  Dataset=paste0(unlist(strsplit(sample_name,"_"))[1],"_",unlist(strsplit(sample_name,"_"))[2])
  Countmatrix=paste0(Scrublet_result_path,sample_name,"barcodes.tsv")
  Countmatrix=read.table(Countmatrix)
  result=cbind(rownames(Countmatrix),sample_scrublet)
  Cell_list=result[which(result[,2]==1),1]
  return(Cell_list)
  
}
Doublets_result=apply(matrix(scrublet_result),1,batch_identify_cellname)
saveRDS(unlist(Doublets_result),file="Datasert_Doublets_cell_list.rds")

Seurat_obj=readRDS(Seurat_obj_path)
Seurat_obj=subset(Seurat_obj,cells=setdiff(rownames(Seurat_obj),Doublets_result))
saveRDS(Seurat_obj,file=paste0(gsub(".rds","",Seurat_obj_path),"_RemoveDoublets.rds"))