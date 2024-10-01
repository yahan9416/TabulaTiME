#' Evaluate the clinical outcomes associated with each cell type in the TCGA (The Cancer Genome Atlas) dataset.
#' @param TCGA_exp_path Path of the merged expression matrix across TCGA cohorts for various cancer types.
#' @param TCGA_clinical_path Path of the clinical information across TCGA cohorts for various cancer types.
#' @param DEG_path Path of DEG file
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
TCGA_exp_path <- args[1]
TCGA_clinical_path <-args[2]
DEG_path <-args[3]


library("survival")
library("survminer")

load(TCGA_exp_path)
load(TCGA_clinical_path)
selecte_cancertype=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC","UVM")
selected_index=match(selecte_cancertype,names(expr))
expr=expr[selected_index]
Survival_info=clin
marker_gene=readRDS(DEG_path)
marker_gene$cluster=as.vector(unlist(marker_gene$cluster))

Survial_cancer_celltype<<-NULL
#batch calculate the coxph for each cell type in every cancer type
cluster_DEgene_effect_one_dataset<-function(clus){
  print(clus)
  temp_marker_gene=marker_gene[which(marker_gene$cluster == clus),]
  DEgenes_cluster=temp_marker_gene$gene[order(temp_marker_gene$avg_log2FC,decreasing = T)]
  
  if(length(DEgenes_cluster)>50){
    DEgenes_cluster=DEgenes_cluster[1:50]
  }else{
    DEgenes_cluster=DEgenes_cluster
  }
  
  batch_get_cor_each_cancertype<-function(temp_cancertype){
    cancer_expMat=t(expr[[temp_cancertype]])
    gsva_score=gsva(cancer_expMat,list(DEgenes_cluster))
    colnames(gsva_score)=substring(colnames(gsva_score),1,12)
    survival_infor=Survival_info
    length(intersect(survival_infor$patient,colnames(gsva_score)))
    sample=intersect(survival_infor$patient,colnames(gsva_score))
    index1=match(sample,survival_infor$patient)
    index2=match(sample,colnames(gsva_score))
    #print(head(survival_infor))
    survival_gsva=cbind(survival_infor[index1,],gsva_score[index2])
    #survival_infor=cbind(survival_infor,gsva_es)
    colnames(survival_gsva)[10]="GSVA_score"
    #cox status: censoring status 1=censored, 2=dead
    survival_gsva$GSVA_score=as.numeric(as.vector(survival_gsva$GSVA_score))
    res.cox <- coxph(Surv(OS, EVENT) ~ GSVA_score, data = survival_gsva)
    temp_result=c(clus,temp_cancertype,summary(res.cox)[[7]],summary(res.cox)[[9]][3])
    
    Survial_cancer_celltype<<-rbind(Survial_cancer_celltype,temp_result)
    
  }
  result=apply(matrix(selecte_cancertype),1,batch_get_cor_each_cancertype)
  
}
result=apply(as.matrix(unique(as.vector(unlist(marker_gene$cluster)))),1,cluster_DEgene_effect_one_dataset)
#result=apply(as.matrix(c("Plasma_JCHAIN","Plasma_CD38","Plasma_IGHM")),1,cluster_DEgene_effect_one_dataset)

colnames(Survial_cancer_celltype)=c("Celltype","Cancertype","coef", "HR","se(coef)", "Z_score", "Pr(>|z|)","P_value")
Survial_cancer_celltype=as.data.frame(Survial_cancer_celltype)
length(unique(Survial_cancer_celltype$Celltype))
Survial_cancer_celltype$P_value=as.numeric(Survial_cancer_celltype$P_value)
Survial_cancer_celltype$Z_score=as.numeric(Survial_cancer_celltype$Z_score)
Survial_cancer_celltype$HR=as.numeric(Survial_cancer_celltype$HR)

##########Visualization############
#heatmap row is cancer type.
survival_heatmap=matrix(Survial_cancer_celltype$Z_score,ncol=length(unique(Survial_cancer_celltype[,1])),byrow=TRUE)
p_value=matrix(Survial_cancer_celltype$P_value,ncol=length(unique(Survial_cancer_celltype[,1])),byrow=TRUE)
colnames(survival_heatmap)=unique(as.vector(unlist(Survial_cancer_celltype$Celltype)))
rownames(survival_heatmap)=unique(Survial_cancer_celltype$Cancertype)
colnames(p_value)=colnames(survival_heatmap)
rownames(p_value)=rownames(survival_heatmap)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
annotation_col=data.frame(Celltype=colnames(survival_heatmap))
rownames(annotation_col)=colnames(survival_heatmap)
annotation_col=as.data.frame(annotation_col)
annotation_row=data.frame(Cancertype=rownames(survival_heatmap))
rownames(annotation_row)=rownames(survival_heatmap)
annotation_row=as.data.frame(annotation_row)
if (!is.null(p_value)){
  ssmt <- p_value< 0.01
  p_value[ssmt] <-'**'
  smt <- p_value >0.01& p_value <0.05
  p_value[smt] <- '*'
  p_value[!ssmt&!smt]<- ''
} else {
  p_value <- F
}


ha_row = rowAnnotation(Treatment=annotation_row$Cancertype)
pdf("TCGA_Cancer_Sep_top50_GSVA_Coxp_Zscore_Pavlue.pdf",height = 10,width = 5.5)
Heatmap(survival_heatmap,c(rep("#414393",3),colorRampPalette(colors = c("#414393","white","#E64536"))(10),rep("#E64536",3)), left_annotation =ha_row,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text( p_value[i, j], x, y, gp = gpar(fontsize = 15, col = "white"))
        })
dev.off()


