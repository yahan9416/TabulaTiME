#' Calculate the signature score for each cell type in order to effectively describe the function of different cell types.
#' In this function, all signarute gene are well known.
#' @param DC_Seu_path         The path of Seurat object only including DC and pDC cells.
#' @param MonoMacro_seu_path  The path of Seurat object only including Monocyte or Macrophage cells.
#' @param output_dir     The path of plot file store

args <- commandArgs(trailingOnly = TRUE)
DC_Seu_path <- args[1]
MonoMacro_seu_path <- args[2]
output_dir <- args[3]

library(MAESTRO)
library(Seurat)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
DC_Seu=readRDS(DC_Seu_path)
DC_signature_score<-function(DC_Seu,output_dir=NULL){
  DefaultAssay(DC_Seu)<-"integrated"
  
  #Signature gene list
  Maturation=c("CD40", "CD80","CD86", "RELB", "CD83")
  Regulatory=c("CD274", "PDCD1LG2", "CD200", "FAS", "ALDH1A2", "SOCS1", "SOCS2")
  TRL_adapotrs=c("MYD88", "MAVS", "TLR9", "TLR8", "TLR7", "TLR6", "TLR5", "TLR4", "TLR3", "TLR2","TLR1")
  Migration=c("CCR7", "MYO1G", "CXCL16", "ADAM8","ICAM1","FSCN1", "MARCKS", "MARCKSL1")
  Th2Response=c("IL4R","IL4I1", "CCL17", "CCL22", "TNFRSF4", "STAT6", "BCL2L1")
  
  meta_infor_name=rownames(DC_Seu@meta.data)
  rownames(DC_Seu@meta.data)=names(Idents(DC_Seu))
  DefaultAssay(DC_Seu)<-"integrated"
  DC_Seu <- AddModuleScore(
    object = DC_Seu,
    features = list(Maturation),
    assay	="RNA",
    ctrl = 5,
    name = 'Maturation'
  )
  DC_Seu <- AddModuleScore(
    object = DC_Seu,
    features = list(Regulatory),
    assay	="RNA",
    ctrl = 5,
    name = 'Regulatory'
  )
  DC_Seu <- AddModuleScore(
    object = DC_Seu,
    features = list(TRL_adapotrs),
    assay	="RNA",
    ctrl = 5,
    name = 'TRL_adapotrs'
  )
  
  DC_Seu <- AddModuleScore(
    object = DC_Seu,
    features = list(Migration),
    assay	="RNA",
    ctrl = 5,
    name = 'Migration'
  )
  DC_Seu <- AddModuleScore(
    object = DC_Seu,
    features = list(Th2Response),
    assay	="RNA",
    ctrl = 5,
    name = 'Th2Response'
  )
  Sign_score= DC_Seu@meta.data
  
  
  Sign_score$source[which(Sign_score$source %in% c("OLK","OSF"))]="Precancerous"
  
  Sign_score_source_celltype<<-NULL
  get_celltype_AverExp<-function(x){
    temp_data=aggregate(Sign_score[,x]~curated_anno_gene+source,Sign_score,mean)
    temp_data$Cell_type=x
    colnames(temp_data)=c("Celltype","source","Signature_score","Signature")
    Sign_score_source_celltype<<-rbind(Sign_score_source_celltype,temp_data)
  }
  result=apply(matrix(c("Maturation1", "Regulatory1", "TRL_adapotrs1", "Migration1", "Th2Response1")),1,get_celltype_AverExp)
  
  #产生行是cell type，列是sig
  cell_type_score<<-NULL
  batch_generate_celltype<-function(x){
    temp_score=Sign_score_source_celltype$Signature_score[which(Sign_score_source_celltype$Celltype ==x)]
    cell_type_score<<-rbind(cell_type_score,temp_score)
    
  }
  result=apply(matrix(unique(DC_Seu@meta.data$curated_anno_gene)),1,batch_generate_celltype)
  colnames(cell_type_score)=paste(rep(c("Maturation1", "Regulatory1", "TRL_adapotrs1", "Migration1", "Th2Response1"),each=3),rep(c("Normal","Precancerous","Tumor"),5),sep = "_")
  rownames(cell_type_score)=unique(DC_Seu@meta.data$curated_anno_gene)
  
  DC_Seu@meta.data=droplevels(DC_Seu@meta.data)
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  pdf("DC_cell_type_source_signateru_col_corr_row.pdf",height = 4,width = 5)
  col_fun = colorRamp2(c( 0,2), c("white", "#276D9F"))
  pheatmap(cell_type_score,cluster_rows = FALSE,cluster_cols=FALSE,scale="row",border_color=NA,gaps_col=seq(3,15,3),color = colorRampPalette(colors = c("#414393","white","#E64536"))(50))
  dev.off()
  
  
}

MonoMacro_signature_score<-function(MonoMacro_seu,output_dir=NULL){
  DefaultAssay(MonoMacro_seu)<-"integrated"
  anti_infla=c("IL1RN","IL10","IL4","IL11","IL13","TGFB1","TNFRSF1A","TNFRSF1B","IL1R2","IL18BP")
  pro_infla=c("IL1B","TNF","CCL2","CCL3","CCL5","CCL7","CCL8","CCL13","CCL17","CCL22")
  M1=c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7" )
  M2=c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")
  Angiogenesis=c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
  Phagocytosis=c("MRC1","CD163","MERTK","C1QB")
  
  
  DefaultAssay(MonoMacro_seu)<-"integrated"
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list("M1"),
    assay	="RNA",
    ctrl = 5,
    name = 'M1_score'
  )
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list(M2),
    assay	="RNA",
    ctrl = 5,
    name = 'M2_score'
  )
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list(Angiogenesis),
    assay	="RNA",
    ctrl = 5,
    name = 'Angiogenesis_score'
  )
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list(Phagocytosis),
    assay	="RNA",
    ctrl = 5,
    name = 'Phagocytosis_score'
  )
  
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list(anti_infla),
    assay	="RNA",
    ctrl = 5,
    name = 'Anti_infla'
  )
  MonoMacro_seu <- AddModuleScore(
    object = MonoMacro_seu,
    features = list(pro_infla),
    assay	="RNA",
    ctrl = 5,
    name = 'Pro_infla'
  )
  
  
  Sign_score= MonoMacro_seu@meta.data
  
  
  Sign_score$source[which(Sign_score$source %in% c("OLK","OSF"))]="Precancerous"
  
  Sign_score_source_celltype<<-NULL
  get_celltype_AverExp<-function(x){
    temp_data=aggregate(Sign_score[,x]~curated_anno_gene+source,Sign_score,mean)
    temp_data$Cell_type=x
    colnames(temp_data)=c("Celltype","source","Signature_score","Signature")
    Sign_score_source_celltype<<-rbind(Sign_score_source_celltype,temp_data)
  }
  result=apply(matrix(c("M1_score", "M2_score", "Angiogenesis_score", "Phagocytosis_score", "Anti_infla","Pro_infla")),1,get_celltype_AverExp)
  
  #产生行是cell type，列是sig
  cell_type_score<<-NULL
  batch_generate_celltype<-function(x){
    temp_score=Sign_score_source_celltype$Signature_score[which(Sign_score_source_celltype$Celltype ==x)]
    cell_type_score<<-rbind(cell_type_score,temp_score)
    
  }
  result=apply(matrix(unique(DC_Seu@meta.data$curated_anno_gene)),1,batch_generate_celltype)
  colnames(cell_type_score)=paste(rep(c("M1_score", "M2_score", "Angiogenesis_score", "Phagocytosis_score", "Anti_infla","Pro_infla"),each=3),rep(c("Normal","Precancerous","Tumor"),6),sep = "_")
  rownames(cell_type_score)=unique(DC_Seu@meta.data$curated_anno_gene)
  
  DC_Seu@meta.data=droplevels(DC_Seu@meta.data)
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  pdf("MonoMacro_cell_type_source_signateru_col_corr_row.pdf",height = 4,width = 5)
  col_fun = colorRamp2(c( 0,2), c("white", "#276D9F"))
  pheatmap(cell_type_score,cluster_rows = FALSE,cluster_cols=FALSE,scale="row",border_color=NA,gaps_col=seq(3,18,3),color = colorRampPalette(colors = c("#414393","white","#E64536"))(50))
  dev.off()
  
  
}
