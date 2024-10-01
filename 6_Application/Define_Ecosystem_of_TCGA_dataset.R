#' Application
#' Define the ecosystem of TCGA samples
#' @param DEG_gene_path Path of merged DEG files across all cell lineags
#' @param TCGA_exp_path Path of the merged expression matrix across TCGA cohorts for various cancer types.
#' @param TCGA_clincal_path Path of the clinical information across TCGA cohorts for various cancer types.
#' @author Ya Han


args <- commandArgs(trailingOnly = TRUE)
DEG_gene_path <- args[1]
TCGA_exp_path <-args[2]
TCGA_clincal_path <-args[3]


#First: Process the DEG of each subtypes
DEG_gene_list=readRDS(DEG_gene_path)
DEG_gene_list=droplevels(DEG_gene_list)
DEG_gene_list=DEG_gene_list[ order(DEG_gene_list$avg_log2FC,decreasing = TRUE),]
Merged_DEG_gene_list<<-list()
process_DEG_list_form<-function(temp_cluster){
  temp_DEG=DEG_gene_list$gene[which(DEG_gene_list$cluster == temp_cluster)]
  if(length(temp_DEG) > 50){
    Merged_DEG_gene_list<<-c(Merged_DEG_gene_list,list(temp_DEG[1:50]))
  }else{
    Merged_DEG_gene_list<<-c(Merged_DEG_gene_list,list(temp_DEG))}
}
result=apply(matrix(unique(DEG_gene_list$cluster)),1,process_DEG_list_form)
names(Merged_DEG_gene_list)=unique(DEG_gene_list$cluster)


########Calculate GSVA score############
library(GSVA)
load(TCGA_exp_path)
load(TCGA_clincal_path)
selecte_cancertype=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC","UVM")
selected_index=match(selecte_cancertype,names(expr))
expr=expr[selected_index]

TCGA_expmat<<-NULL
batch_get_cor_each_cancertype<-function(temp_cancertype){
  cancer_expMat=t(expr[[temp_cancertype]])  
  TCGA_expmat<<-cbind(TCGA_expmat,cancer_expMat)
}
result=apply(matrix(selecte_cancertype),1,batch_get_cor_each_cancertype)

TCGA_comb_gsva=gsva(TCGA_expmat,DEG_gene_list)

########Cluster and visualization###########
pdf("TCGA_LineDEG_heatmap_cluster.pdf",height = 7,width = 18)
pheatmap(TCGA_comb_gsva,cluster_cols=TRUE,cluster_rows=TRUE,scale="none",color =c(colorRampPalette(colors = c("#460054","#399982","#FBE721"))(20)),border_color=NA,clustering_method="average",cutree_cols=3,show_colnames=FALSE)
dev.off()  

