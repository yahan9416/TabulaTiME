#' Application: Benchmarking the performance of cell type annotations in scRNA-seq dataset between signatures identified by TabulaTIME and Cellstates
#' @param Seuratobj_path Path of test datasets
#' @param TabulaTIME_DEG_path Path of merged DEG file of TabulaTIME
#' @param TabulaTIME_DEG_SCINA_path Path of Tabula_DEG_SCINA_annotation_signature.csv
#' @param Cell_69cellstates_path Path of Cell_Altas_TCGAsample_69cellstates.txt
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Seuratobj_path <- args[1]
TabulaTIME_DEG_path <-args[2]
TabulaTIME_DEG_SCINA_path <-args[3]
Cell_69cellstates_path<-args[4]


Seuratobj=readRDS(Seuratobj_path)
Tabula_DEG=readRDS(TabulaTIME_DEG_path)
#######Cell state signature score in addation datasets###########
batch_process_lineage_score<-function(temp_lineage){
  temp_cells=rownames(Seuratobj@meta.data)[which(Seuratobj@meta.data$assign.curated_major == temp_lineage)]
  lineage_seu=subset(Seuratobj,cells=temp_cells)
  lineage_seu <- NormalizeData(lineage_seu)
  lineage_seu <- FindVariableFeatures(lineage_seu, selection.method = "vst", nfeatures = 2000)
  lineage_seu <- ScaleData(lineage_seu)
  lineage_seu <- RunPCA(lineage_seu)
  lineage_seu <- RunUMAP(lineage_seu, dims = 1:10)
  lineage_seu <- FindNeighbors(lineage_seu, dims = 1:10)
  lineage_seu <- FindClusters(lineage_seu, resolution = 0.5)
  
  TabulaTime_sign=Tabula_DEG[which(TabulaTiME_lineag_name == temp_lineage)]

  
  library(GSVA)
  expmat=GetAssayData(lineage_seu)
  expmat=as.matrix(expmat)
  Signature_score_TabulaTIME=gsva(expmat,TabulaTime_sign,method="ssgsea")
  Signature_score_cell69state=gsva(expmat,Cell_69cellstates_list[1:5],method="ssgsea")
  
  lineage_seu@meta.data$seurat_clusters=as.vector(unlist(lineage_seu@meta.data$seurat_clusters))
  cell_type_score<<-NULL
  batch_generate_celltype<-function(x){
    temp_score=rownames(lineage_seu@meta.data)[which(lineage_seu@meta.data$seurat_clusters ==x)]
    celltype_sig_TabulaTIME=rowMeans(Signature_score_TabulaTIME[,match(temp_score,colnames(Signature_score_TabulaTIME))])
    celltype_sig_Cellstate=rowMeans(Signature_score_cell69state[,match(temp_score,colnames(Signature_score_cell69state))])
    
    cell_type_score<<-rbind(cell_type_score,celltype_sig_TabulaTIME,celltype_sig_Cellstate)
    
  }
  result=apply(matrix(unique(lineage_seu@meta.data$seurat_clusters)),1,batch_generate_celltype)
  
  
  rownames(cell_type_score)=unique(lineage_seu@meta.data$seurat_clusters)
}
result=apply(matrix(unique(TabulaTiME_lineag_name)),1,batch_process_lineage_score)



######Accuracy of cell type prediction###############
#####TabulaTiME######
Tabula_DEG=readRDS(Tabula_DEG_path)
Tabula_DEG_matrix<<-NULL
result=lapply(Tabula_DEG,function(x){Tabula_DEG_matrix<<-cbind(Tabula_DEG_matrix,unlist(x))})
colnames(Tabula_DEG_matrix)=names(Tabula_DEG)

library('SCINA')
signatures=preprocess.signatures('TabulaTIME_DEG_SCINA_path')

##########Cell_69cellstates#############
Cell_69cellstates=read.table(Cell_69cellstates_path,header = TRUE,sep="\t")
Cell_69cellstates[8,2]="SEPT1"
Cell_69cellstates$Celltype_state=paste0(unlist(lapply(strsplit(Cell_69cellstates[,1]," "),function(x) x[1])),"_",Cell_69cellstates[,3])
Cell_69cellstates$Celltype_state=gsub("Monocytes/Macrophages","MonoMacro",Cell_69cellstates$Celltype_state)
#Transfer into list formation
Cell_69cellstates_list<-NULL
batch_get_DEG_list<-function(clus){
  temp_marker_gene=Cell_69cellstates$Gene[which(Cell_69cellstates$Celltype_state == clus)]
  Cell_69cellstates_list<<-c(Cell_69cellstates_list,list(temp_marker_gene))
}
result=apply(matrix(unique(Cell_69cellstates$Celltype_state)),1,batch_get_DEG_list)
names(Cell_69cellstates_list)=unique(Cell_69cellstates$Celltype_state)

###########SCINA Predict##################
library(MAESTRO)
library(Seurat)

expmat=GetAssayData(Seuratobj) #
expmat=as.matrix(expmat)

results = SCINA(expmat, signatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')

Cell69_results = SCINA(expmat, Cell_69cellstates_list, max_iter = 100, convergence_n = 10, 
                       convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')

predict_result=data.frame(CellID=colnames(expmat),Pre_celltype=results[["cell_labels"]])
predict_result$Orignal_celltype=Seuratobj@meta.data$assign.curated
predict_result$Cell69_PreCelltype=Cell69_results[["cell_labels"]]

Accuracy_TabulaTiME=(length(which(predict_result$Orignal_celltype == predict_result$TabulaTiME_PreCelltype_major))+length(which(predict_result$TabulaTiME_PreCelltype_major == "unknown")))/length(predict_result$Orignal_celltype)

#temp_predict_result=predict_result[-1*which(predict_result$Cell69_PreCelltype_major == "unknown"),]
Accuracy_Cell69=(length(which(predict_result$Cell69_PreCelltype_major == predict_result$Orignal_celltype))+length(which(predict_result$Cell69_PreCelltype_major == "unknown")))/length(predict_result$Orignal_celltype)
