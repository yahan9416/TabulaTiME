library(Seurat)
library(ggplot2)
library(MAESTRO)


#The first things is to find which dataset with source information and need to Down-sampling 
Tumor_file=list.files("./Seurat_obj",recursive=T,full.names = T)
batch_calculate_TPM<-function(x){
  seurat_obj=readRDS("./Data/GEO/NSCLC_GSE162498/NSCLC_GSE162498_res.rds")
  count=GetAssayData(seurat_obj$RNA,slot ="counts")
  cal_TPM <- RNACountToTPM(count, idType = "SYMBOL", organism = "GRCh38")
  log_TPM=log2((cal_TPM/10)+1)
  x=gsub("_res","_TPM_exp",x)
  setwd("./scRNAseq_1OX/TPM")
  saveRDS(log_TPM,file=x)
  
}
result=apply(as.matrix(Tumor_file),1,batch_calculate_TPM)
