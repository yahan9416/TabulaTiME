#Application
#define the ecosystem of TCGA samples

#First: Process the DEG of each subtypes
DEG_gene_list=readRDS("./Endothelial/Basic_analysis/Endo_celltype_DEG.rds")
DEG_gene_list=droplevels(DEG_gene_list)
DEG_gene_list=DEG_gene_list[ order(DEG_gene_list$avg_log2FC,decreasing = TRUE),]
Endo_DEG_gene_list<<-list()
process_DEG_list_form<-function(temp_cluster){
  temp_DEG=DEG_gene_list$gene[which(DEG_gene_list$cluster == temp_cluster)]
  if(length(temp_DEG) > 50){
    Endo_DEG_gene_list<<-c(Endo_DEG_gene_list,list(temp_DEG[1:50]))
  }else{
    Endo_DEG_gene_list<<-c(Endo_DEG_gene_list,list(temp_DEG))}
}
result=apply(matrix(unique(DEG_gene_list$cluster)),1,process_DEG_list_form)
names(Endo_DEG_gene_list)=unique(DEG_gene_list$cluster)
##########Fibroblast###
marker_gene=readRDS("./Fibroblast/Basic_analysis/Fibroblast_celltype_DEG.rds")
marker_gene=marker_gene[which(marker_gene$avg_log2FC > 0),]
marker_gene$cluster=as.vector(unlist(marker_gene$cluster))

Fibro_DEG_gene_list<<-list()
process_DEG_list_form<-function(temp_cluster){
  temp_DEG=marker_gene$gene[which(marker_gene$cluster == temp_cluster)]
  if(length(temp_DEG) > 50){
    Fibro_DEG_gene_list<<-c(Fibro_DEG_gene_list,list(temp_DEG[1:50]))
  }else{
    Fibro_DEG_gene_list<<-c(Fibro_DEG_gene_list,list(temp_DEG))}
}
result=apply(matrix(unique(marker_gene$cluster)),1,process_DEG_list_form)
names(Fibro_DEG_gene_list)=unique(marker_gene$cluster)

DEG_gene_list=c(Endo_DEG_gene_list,Fibro_DEG_gene_list)
cell_type_lineage=rbind(cbind("Endo",names(Endo_DEG_gene_list)),cbind("Fibro",names(Fibro_DEG_gene_list)))

########Calculate GSVA score############
library(GSVA)
load('/fs/home/hanya/Project/Survival_cohorts/TCGAexpr.RData')
load("/fs/home/hanya/Project/Survival_cohorts/TCGAclin.RData")
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
saveRDS(TCGA_comb_gsva,file="./Apply_bulk/All_lineage_DEG_together_TCGA.rds")

########Cluster and visualization###########
pdf("TCGA_LineDEG_heatmap_cluster.pdf",height = 7,width = 18)
pheatmap(TCGA_comb_gsva,cluster_cols=TRUE,cluster_rows=TRUE,scale="none",color =c(colorRampPalette(colors = c("#460054","#399982","#FBE721"))(20)),border_color=NA,clustering_method="average",cutree_cols=3,show_colnames=FALSE)
dev.off()  

