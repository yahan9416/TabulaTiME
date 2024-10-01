#' Inferring the probability of metastasis
#' @param Plus_preM_gene_path Path of PLUS_TCGA_191_metastasis_predicted_genes.csv
#' @param Datainfo_path Information of scRNA-seq datasets 
#' @author Ya Han


args <- commandArgs(trailingOnly = TRUE)
Plus_preM_gene_path <- args[1]
TCGA_exp_path <-args[2]

#Load the Pus predefined the metastasis gene list
Plus_preM_gene=as.vector(unlist(read.csv(Plus_preM_gene_path,header = FALSE)))
Plus_preM_gene=list(Plus_preM_gene)
#
library(GSVA)
load(TCGA_exp_path)
selecte_cancertype=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC","UVM")
selected_index=match(selecte_cancertype,names(expr))
expr=expr[selected_index]



TCGA_expmat<<-NULL
batch_get_cor_each_cancertype<-function(temp_cancertype){
  cancer_expMat=t(expr[[temp_cancertype]])  
  TCGA_expmat<<-cbind(TCGA_expmat,cancer_expMat)
}
result=apply(matrix(selecte_cancertype),1,batch_get_cor_each_cancertype)

TCGA_comb_gsva=gsva(TCGA_expmat,Plus_preM_gene)
