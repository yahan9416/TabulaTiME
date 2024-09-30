#' Generating the TPM from Seurat object
#' @param File_path Path of Seurat_obj
#' @author Ya Han

library(Seurat)
library(ggplot2)
library(MAESTRO)

#batch generate the TPM 
 args <- commandArgs(trailingOnly = TRUE) #The direction or file_path of Seurat_obj
 File_path <- args[1]
 
#The first things is to find which dataset with source information and need to Down-sampling 
Tumor_file=list.files(File_path,recursive=T,full.names = T)
Tumor_file=Tumor_file[grep("Seurat_obj",Tumor_file)]
batch_calculate_TPM<-function(x){
  seurat_obj=readRDS(x)
  count=GetAssayData(seurat_obj,slot ="counts")
  cal_TPM <- RNACountToTPM(count, idType = "SYMBOL", organism = "GRCh38")
  log_TPM=log2((cal_TPM/10)+1)
  x=gsub("_Seurat_obj","_TPM_exp",x) #Convert the Seurat_obj name to the TPM file name
  saveRDS(log_TPM,file=x)
}
result=apply(as.matrix(Tumor_file),1,batch_calculate_TPM)
