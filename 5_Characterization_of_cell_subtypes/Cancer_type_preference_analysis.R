Datainfo=read.table("./TISCH_V1V2_10X_data_patient_cell.txt",header = TRUE,sep="\t")
Treatment_naive=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "None")]
Immthe_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Immunotherapy")]
Chemo_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Chemotherapy")]

Cancer_type_class=read.table("./TISCH_cancertype_classify.txt",header = TRUE,sep="\t")


meta.data$Sample=as.vector(unlist(lapply(strsplit(rownames(meta.data),"\\|"),function(x)x[1])))

cell_index=c(which(meta.data$Cancer_type %in% c("Breast","Live","Oral","Kidney","Ovary","PBMC")),intersect(which(meta.data$Cancer_type %in% setdiff(unique(meta.data$Cancer_type),c("Breast","Live","Oral","Kidney","Ovary","PBMC"))),which(meta.data$Source %in% c("Tumor"))))
cell_index=unique(c(which(meta.data$Cancer_type %in% c("Breast","Live","Oral","Kidney","Ovary","PBMC")),which(meta.data$Source %in% c("Tumor"))))

meta.data=meta.data[cell_index,]
stat_result=table(meta.data$Dataset,meta.data$curated_anno)
Propor_ICB_data<<-NULL
Convert_Dataset_Cancertype<-function(temp_cancer){
  mini_cancer=Cancer_type_class[which(Cancer_type_class$Major_Cancer == temp_cancer),1]
  Dataset=unique(meta.data$Dataset[which(meta.data$Cancer_type %in% mini_cancer)])
  ICB_data=intersect(Immthe_Dataset,Dataset)
  Non_ICB_data=setdiff(Dataset,ICB_data)
  
  if(length(ICB_data) >1){
    temp_icb=c(temp_cancer,"ICB_Treatment",length(which(meta.data$Dataset %in% ICB_data)),colSums(stat_result[match(ICB_data,rownames(stat_result)),])/sum(stat_result[match(ICB_data,rownames(stat_result)),]))
  }
  if(length(ICB_data) == 1){
    temp_icb=c(temp_cancer,"ICB_Treatment",length(which(meta.data$Dataset %in% ICB_data)),stat_result[match(ICB_data,rownames(stat_result)),]/sum(stat_result[match(ICB_data,rownames(stat_result)),]))
  }
  if(length(Non_ICB_data) >1){
    temp_nonicb=c(temp_cancer,"Treatment_naive",length(which(meta.data$Dataset %in% Non_ICB_data)),colSums(stat_result[match(Non_ICB_data,rownames(stat_result)),])/sum(stat_result[match(Non_ICB_data,rownames(stat_result)),]))
  }
  if(length(Non_ICB_data) == 1){
    temp_nonicb=c(temp_cancer,"Treatment_naive",length(which(meta.data$Dataset %in% Non_ICB_data)),stat_result[match(Non_ICB_data,rownames(stat_result)),]/sum(stat_result[match(Non_ICB_data,rownames(stat_result)),]))
  }
  
  if(length(Non_ICB_data) >= 1 && length(ICB_data)>=1){
    Propor_ICB_data<<-rbind(Propor_ICB_data,rbind(temp_icb,temp_nonicb))}
  if(length(Non_ICB_data) >= 1 && length(ICB_data)==0){
    Propor_ICB_data<<-rbind(Propor_ICB_data,temp_nonicb)}
  if(length(Non_ICB_data) == 0 && length(ICB_data)>=1 ){
    Propor_ICB_data<<-rbind(Propor_ICB_data,temp_icb)}
  
}
result=apply(matrix(unique(Cancer_type_class$Major_Cancer)),1,Convert_Dataset_Cancertype)

rownames(Propor_ICB_data)=paste(Propor_ICB_data[,1],Propor_ICB_data[,2],sep="|")

library(pheatmap)
annotation_row=cbind(Cancer_type_class[match(Propor_ICB_data[,1],Cancer_type_class[,2]),5],Propor_ICB_data[,2])
annotation_row[which(annotation_row[,1] == "Healthy_derived" ),2]="Healthy_derived"
annotation_row=as.data.frame(annotation_row)
rownames(annotation_row)=rownames(Propor_ICB_data)
#Cell_number_cancer_ICB=data.frame(Cell_number=as.numeric(Present_ICB_data[,3]))
Cell_number_cancer_ICB=data.frame(Cell_number=log2(as.numeric(Propor_ICB_data[,3])))
rownames(Cell_number_cancer_ICB)=rownames(Propor_ICB_data)

celltype_order=colnames(Propor_ICB_data)[-1*c(1,2,3)]
Propor_ICB_data=matrix(as.numeric(unlist(Propor_ICB_data[,-1*c(1,2,3)])),byrow=FALSE,nrow=dim(Propor_ICB_data)[1])
colnames(Propor_ICB_data)=celltype_order
rownames(Propor_ICB_data)=rownames(annotation_row)



library(ComplexHeatmap)
ha_row = rowAnnotation(ICB=annotation_row$V2,
                       Originated=annotation_row$V1,
                       Log_Cell_number = anno_barplot(Cell_number_cancer_ICB$Cell_number),
                       annotation_name_rot = 45,col=list(ICB=c("Treatment_naive"="#91E7FC","Healthy_derived"="#B3FC9A","ICB_Treatment"="#D1B9DD"),Originated=c("Epithelial_Originated"="#91E7FC","Healthy_derived"="#F790E7","Neuroendocrine"="#716EE0", "Sarcomas"="#88F3A9","Melanoma"="#9BD210"))
)



#match(rownames(Propor_ICB_data),paste0(WithinCancertype$Cancertype,"|",WithinCancertype$Treatment))
p_value=matrix(rep(1,length(Propor_ICB_data)),ncol=ncol(Propor_ICB_data))
colnames(p_value)=colnames(Propor_ICB_data)
rownames(p_value)=rownames(Propor_ICB_data)
p_value[,7]=WithinCancertype$FDR[match(rownames(Propor_ICB_data),paste0(WithinCancertype$Cancertype,"|",WithinCancertype$Treatment))]
if (!is.null(p_value)){
  ssmt <- p_value< 0.01
  p_value[ssmt] <-'**'
  smt <- p_value >0.01& p_value <0.05
  p_value[smt] <- '*'
  p_value[!ssmt&!smt]<- ''
} else {
  p_value <- F
}

pdf("Fibro_CancerType_Dataset_Proportion_Cell_proportion_addWithinCancefr_Sign.pdf",height = 10,width = 8)
Heatmap(Propor_ICB_data,c(colorRampPalette(colors = c("white","#E64536"))(20),rep("#E64536",30)), left_annotation =ha_row,row_split = data.frame(annotation_row$V2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text( p_value[i, j], x, y, gp = gpar(fontsize = 17, col = "grey"))
        })
dev.off()
