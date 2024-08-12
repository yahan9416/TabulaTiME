library(MAESTRO)
library(Seurat)
library(ggplot2)
setwd("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/P13/Integrated/")
Oral_ST_P13_Int=readRDS("Oral_ST_P13_Integrated_Seurat.rds")
# plot the expression of gene pair
DefaultAssay(Oral_ST_P13_Int)<-"SCT"
exp=GetAssayData(Oral_ST_P13_Int)
index_ccl5=match("CCL5",rownames(exp))
index_ccr5=match("CCR5",rownames(exp))

CCL5_CCR5=as.numeric(exp[index_ccl5,])*as.numeric(exp[index_ccr5,])
Oral_ST_P13_Int@meta.data$CCL5_CCR5=CCL5_CCR5

Oral_ST_P13_Int[["CCL5_CCR5"]] <- CCL5_CCR5

pdf("Oral_ST_P13_Int_SCT_CCL5_CCR5.pdf",width = 12,height = 8)
SpatialFeaturePlot(Oral_ST_P13_Int, features = "CCL5_CCR5")
dev.off()
#
#############calcaluate interaction including the neighborhood cell effect########
P13_NorOLK=readRDS("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/P13/Oral_ST_P13Normal_seurat.rds")
P13_Tumor=readRDS("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/P13/Oral_ST_P13Tumor_seurat.rds")
loc_NorOLK=GetTissueCoordinates(P13_NorOLK, scale = NULL)
loc_Tumor=GetTissueCoordinates(P13_Tumor, scale = NULL)
loc_Tumor[,1]=as.numeric(loc_Tumor[,1])
loc_Tumor[,2]=as.numeric(loc_Tumor[,2])
loc_NorOLK[,1]=as.numeric(loc_NorOLK[,1])
loc_NorOLK[,2]=as.numeric(loc_NorOLK[,2])
rownames(loc_NorOLK)<-paste(rownames(loc_NorOLK),"_1",sep="")
rownames(loc_Tumor)<-paste(rownames(loc_Tumor),"_2",sep="")

#test cell  rownames(Oral_ST_P13_Int[["CCL5_CCR5"]])[which(Oral_ST_P13_Int[["CCL5_CCR5"]] >2)]

index=order(Oral_ST_P13_Int@meta.data$CCL5_CCR5,decreasing = TRUE)[1:5]
Oral_ST_P13_Int@meta.data$CCL5_CCR5[index]
as.numeric(loc_NorOLK[match("TAATAGAACAGAGTTA-1_1",rownames(loc_NorOLK)),])
as.numeric(loc_NorOLK[match("AGCGACAGGAACGGTC-1_1",rownames(loc_NorOLK)),])

calculate_cell_interaction<-function(temp_cell){
  source_index=strsplit(temp_cell,"_")[[1]][2]
  if(source_index == "1"){
    cell_loca=as.numeric(loc_NorOLK[match(temp_cell,rownames(loc_NorOLK)),])
    cell_neig_row=intersect(which(loc_NorOLK[,1] <= (cell_loca[1]+200)),which(loc_NorOLK[,1] >= (cell_loca[1]-200)))
    cell_neig_col=intersect(which(loc_NorOLK[,2] <= (cell_loca[2]+200)),which(loc_NorOLK[,2] >= (cell_loca[2]-200)))
    index=intersect(cell_neig_col,cell_neig_row)
    print(length(index))
    neig_index=match(setdiff(rownames(loc_NorOLK)[index],temp_cell),colnames(exp))
    temp_cell_index=match(temp_cell,colnames(exp))
    CCL5_CCR5_tempcell=as.numeric(exp[index_ccl5,temp_cell_index])*as.numeric(exp[index_ccr5,temp_cell_index])+0.5*sum(as.numeric(exp[index_ccl5,neig_index])*as.numeric(exp[index_ccr5,neig_index]))
    return(CCL5_CCR5_tempcell)
  }
  if(source_index == "2"){
    cell_loca=as.numeric(loc_Tumor[match(temp_cell,rownames(loc_Tumor)),])
    cell_neig_row=intersect(which(loc_Tumor[,1] <= (cell_loca[1]+200)),which(loc_Tumor[,1] >= (cell_loca[1]-200)))
    cell_neig_col=intersect(which(loc_Tumor[,2] <= (cell_loca[2]+200)),which(loc_Tumor[,2] >= (cell_loca[2]-200)))
    index=intersect(cell_neig_col,cell_neig_row)
    print(length(index))
    neig_index=match(setdiff(rownames(loc_Tumor)[index],temp_cell),colnames(exp))
    temp_cell_index=match(temp_cell,colnames(exp))
    CCL5_CCR5_tempcell=as.numeric(exp[index_ccl5,temp_cell_index])*as.numeric(exp[index_ccr5,temp_cell_index])+0.5*sum(as.numeric(exp[index_ccl5,neig_index])*as.numeric(exp[index_ccr5,neig_index]))
    return(CCL5_CCR5_tempcell)
  }
  
}
result=apply(matrix(rownames(Oral_ST_P13_Int[["CCL5_CCR5"]])),1,calculate_cell_interaction)
Oral_ST_P13_Int[["CCL5_CCR5"]]=as.numeric(unlist(result))

pdf("Oral_ST_P13_Int_SCT_CCL5_CCR5_Neig6.pdf",width = 12,height = 8)
SpatialFeaturePlot(Oral_ST_P13_Int, features = "CCL5_CCR5")
dev.off()
