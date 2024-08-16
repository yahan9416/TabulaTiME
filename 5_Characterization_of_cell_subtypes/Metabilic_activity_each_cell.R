#Metabolic


#First: apply magic to imputation

library(Rmagic)
library(phateR)
library(scater)
library(MAESTRO)
library(Seurat)
library(ggplot2)

Seurat_obj=readRDS(temp_file)
DefaultAssay(Seurat_obj$RNA)<-"RNA"
count=GetAssayData(Seurat_obj$RNA,slot ="counts")
cal_TPM <- RNACountToTPM(count, idType = "SYMBOL", organism = "GRCh38")
log_TPM=log2((cal_TPM/10)+1)
log_TPM=as.matrix(log_TPM)

if(ncol(log_TPM)>5000){knn=15}
  data_MAGIC <- magic(t(log_TPM), knn=knn,n.jobs=10,seed=1, genes="all_genes" )
  imputation_TPM=data_MAGIC$result
  imputation_TPM=as.matrix(imputation_TPM)
  imputation_TPM=t(imputation_TPM)

  
######
 
  load("./Metabolic_result/KEGG_Module/Module_gene_process/KEGG_Module_genelist_nameinfo.RData")
  load("./Metabolic_result/MetaCyc/MetaCyc_Homo_onlyKeepHS_pathway_gene.RData")
  load("./Metabolic_result/Enzyme_relationship_corr_uni.RData")
  Enzyme_or_matrix=imputation_TPM[na.omit(match(unlist(Enzyme_or_list_uni),rownames(imputation_TPM))),]
  gene_name_or<<-NULL
  Enzyme_or_matrix_correct<<-NULL
  batch_process_Enzyme_or_gene_exp<-function(x){
    index=na.omit(match(x,rownames(Enzyme_or_matrix)))
    if(length(index) == 1){
      gene_name_or<<-c(gene_name_or,rownames(Enzyme_or_matrix)[index])
      Enzyme_or_matrix_correct<<-rbind(Enzyme_or_matrix_correct,Enzyme_or_matrix[index,])
    }
    if(length(index) > 1){
      temp_expr=Enzyme_or_matrix[na.omit(match(x,rownames(Enzyme_or_matrix))),]
      gene_name_or<<-c(gene_name_or,rownames(Enzyme_or_matrix)[index])
      result=apply(temp_expr,2,function(y) max(y))
      result=matrix(rep(as.vector(unlist(result)),length(index)),byrow=TRUE,nrow=length(index))
      Enzyme_or_matrix_correct<<-rbind(Enzyme_or_matrix_correct,result)
    }
  }
  result=lapply(Enzyme_or_list_uni,batch_process_Enzyme_or_gene_exp)
  rownames(Enzyme_or_matrix_correct)<-gene_name_or
  
################################################
  #Correct and enzyme
  Enzyme_and_matrix=imputation_TPM[na.omit(match(unlist(Enzyme_and_list_uni),rownames(imputation_TPM))),]
  gene_name_and<<-NULL
  Enzyme_and_matrix_correct<<-NULL
  batch_process_Enzyme_or_gene_exp<-function(x){
    index=na.omit(match(x,rownames(Enzyme_and_matrix)))
    if(length(index) == 1){
      gene_name_and<<-c(gene_name_and,rownames(Enzyme_and_matrix)[index])
      Enzyme_and_matrix_correct<<-rbind(Enzyme_and_matrix_correct,Enzyme_and_matrix[index,])
      #return(Enzyme_or_matrix[index,])
    }
    if(length(index) > 1){
      temp_expr=Enzyme_and_matrix[na.omit(match(x,rownames(Enzyme_and_matrix))),]
      gene_name_and<<-c(gene_name_and,rownames(Enzyme_and_matrix)[index])
      result=apply(temp_expr,2,function(y) max(y))
      result=matrix(rep(as.vector(unlist(result)),length(index)),byrow=TRUE,nrow=length(index))
      Enzyme_and_matrix_correct<<-rbind(Enzyme_and_matrix_correct,result)
      #return(result)
    }
  }
  result=lapply(Enzyme_and_list_uni,batch_process_Enzyme_or_gene_exp)
  rownames(Enzyme_and_matrix_correct)<-gene_name_and
  #其实pathway也不是必须要修改，因为没有的大家都咩有，所有也就无所谓了
  
  replace_gene_index=unique(c(na.omit(match(unlist(Enzyme_and_list_uni),rownames(imputation_TPM))),na.omit(match(unlist(Enzyme_or_list_uni),rownames(imputation_TPM)))))
  
  imputation_TPM=imputation_TPM[-1*replace_gene_index,]
  imputation_TPM=rbind(imputation_TPM,Enzyme_and_matrix_correct)
  imputation_TPM=rbind(imputation_TPM,Enzyme_or_matrix_correct)
  
  pathway_activity<-gsva(imputation_TPM,Module_gene_list,parallel.sz=20)
  saveRDS(pathway_activity,file=paste0("KEGG_metabolic_Module_addEnzy_GSVA_activity.rds"))
  