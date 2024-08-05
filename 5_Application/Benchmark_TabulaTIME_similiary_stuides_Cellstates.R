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
Tabula_DEG=readRDS("./Lineage_RNAassay_DEG_top50.rds")
Tabula_DEG_matrix<<-NULL
result=lapply(Tabula_DEG,function(x){Tabula_DEG_matrix<<-cbind(Tabula_DEG_matrix,unlist(x))})
colnames(Tabula_DEG_matrix)=names(Tabula_DEG)
write.table(Tabula_DEG_matrix,file="Tabula_DEG_SCINA_annotation_signature.csv",col.names = TRUE,row.names = FALSE,sep=",",quote = FALSE)

library('SCINA')
signatures=preprocess.signatures('./Tabula_DEG_SCINA_annotation_signature.csv')

##########Cell_69cellstates#############
Cell_69cellstates=read.table("/fs/home/hanya/Project/TME_Immune_difference/NatCancer_revision/Apply_bulk/Benchmark_othermethods/Cell_TenEcotypes/Cell_Altas_TCGAsample_69cellstates.txt",header = TRUE,sep="\t")
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
Seuratobj=readRDS(paste0("./BRCA_GSE176078/",file_list[grep("_res.rds",file_list)]))
Seuratobj=Seuratobj$RNA
expmat=GetAssayData(Seuratobj) #这里面的值已经是normalization过后的
expmat=as.matrix(expmat)
#expmat_raw=Read10X_h5("/fs/home/wangyuting/biosoft/django/wyt_test/TISCH1/static/data/OV_GSE147082/OV_GSE147082_TPM.h5")
#expmat_raw[1:5,match(c("GSM4416534@ACCCTCTCGCTA", "GSM4416534@ATTAGGAGACTG"),colnames(expmat_raw))]



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
