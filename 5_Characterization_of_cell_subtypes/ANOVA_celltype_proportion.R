#' Evaluting the source distribution of cell types within each lineage by ANOVA test 
#' @param meta_info_path Path of meta.data
#' @param Datainfo_path Information of scRNA-seq datasets 
#' @author Ya Han


args <- commandArgs(trailingOnly = TRUE)
meta_info_path <- args[1]
Datainfo_path <-args[2]

meta_info=readRDS(meta_info_path)
meta_info$Sample=as.vector(unlist(lapply(strsplit(rownames(meta_info),"\\|"),function(x)x[1])))
#meta_info$Sample=as.vector(unlist(lapply(strsplit(meta_info$Sample,"\\@"),function(x)x[2])))

#Proportion
cluster_sample=table(meta_info$Sample,meta_info$curated_anno)
source_pro=t(apply(cluster_sample,1,function(x) x/sum(x)))
sample_pro=data.frame(Sample=rep(rownames(source_pro),dim(source_pro)[2]),Proportion=as.numeric(as.vector(unlist(source_pro))),Celltype=rep(colnames(source_pro),each=dim(source_pro)[1]))
sample_pro$Source=meta_info$Source[match(sample_pro$Sample,meta_info$Sample)]
sample_pro$Dataset=meta_info$Dataset[match(sample_pro$Sample,meta_info$Sample)]

#Dataset info
Datainfo=read.table(Datainfo_path,header = TRUE,sep="\t")
Treatment_naive=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "None")]
Immthe_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Immunotherapy")]
Chemo_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Chemotherapy")]

library(ggplot2)
library(ggpubr)
Without_Treatment=sample_pro[which(sample_pro$Dataset %in% Treatment_naive),]
Without_Treatment$Source=factor(Without_Treatment$Source,levels=c("Blood","Normal","Precancerous","Tumor","Metastatic" ))
Without_Treatment=Without_Treatment[which(!is.na(Without_Treatment$Source)),]


batch_celltype_source_distribution<-function(temp_celltype){
  temp_data=Without_Treatment[which(Without_Treatment$Celltype == temp_celltype ),]
  aov_result=aov(Proportion~Source,data=temp_data)
  P_value=unlist(summary(aov_result)[1])[9]
  return(c(temp_celltype,P_value))
}
result=apply(matrix(unique(Without_Treatment$Celltype)),1,batch_celltype_source_distribution)
result=t(result)
aov_res=as.data.frame(result)
colnames(aov_res)=c("Celltype","P_value")
aov_res$P_value=as.numeric(aov_res$P_value)
aov_res$FDR=p.adjust(aov_res$P_value,method="BH")

