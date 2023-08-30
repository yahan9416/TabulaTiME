###Calculate the probability of metastasis

#Load the Pus predefined the metastasis gene list
Plus_preM_gene=as.vector(unlist(read.csv("./PLUS_TCGA_191_metastasis_predicted_genes.csv",header = FALSE)))
Plus_preM_gene=list(Plus_preM_gene)
#我们可不可以尝试对一个slide 来说对malignant cell 计算metastasis score 以及对fibroblast cell 计算CTHRC1 score
library(GSVA)
load('/fs/home/hanya/Project/Survival_cohorts/TCGAexpr.RData')
selecte_cancertype=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC","UVM")
selected_index=match(selecte_cancertype,names(expr))
expr=expr[selected_index]


#在测试后发现在一个癌型的expression matrix 中和两个cancer type的一起计算时，得到的结果是不一样的，所以把所有的放在一起计算

TCGA_expmat<<-NULL
batch_get_cor_each_cancertype<-function(temp_cancertype){
  cancer_expMat=t(expr[[temp_cancertype]])  
  TCGA_expmat<<-cbind(TCGA_expmat,cancer_expMat)
}
result=apply(matrix(selecte_cancertype),1,batch_get_cor_each_cancertype)

TCGA_comb_gsva=gsva(TCGA_expmat,Plus_preM_gene)
