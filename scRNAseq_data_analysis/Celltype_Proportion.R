#Estimate the proportion of each cell type
SeuratObj@meta.data=meta_info
meta_info=Seurat_obj@meta.data
meta_info$curated_anno=droplevels(meta_info$curated_anno)
meta_info$Sample=as.vector(unlist(lapply(strsplit(rownames(meta_info),"\\|"),function(x)x[1])))

#Divided dataset into treatment naive and therapy dataset
Datainfo=read.table("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/TISCH_V1V2_10X_data_patient_cell1.txt",header = TRUE,sep="\t")
Treatment_naive=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "None")]
Immthe_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Immunotherapy")]
Chemo_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Chemotherapy")]

cluster_sample=table(meta_info$Sample,meta_info$curated_anno)
source_pro=t(apply(cluster_sample,1,function(x) x/sum(x)))
sample_pro=data.frame(Sample=rep(rownames(source_pro),dim(source_pro)[2]),Proportion=as.vector(source_pro),Cell_type=rep(colnames(source_pro),each=dim(source_pro)[1]))

sample_pro$Tissue=meta_info$Tissue[match(sample_pro$Sample,meta_info$Sample)]
sample_pro$Source=meta_info$Source[match(sample_pro$Sample,meta_info$Sample)]
sample_pro$Dataset=meta_info$Dataset[match(sample_pro$Sample,meta_info$Sample)]



#######Cell type proportion across Teatment naive dataset#####
Without_Treatment=sample_pro[which(sample_pro$Dataset %in% Treatment_naive),]
Without_Treatment$Source=factor(Without_Treatment$Source,levels=c("Blood","Normal","Precancerous","Tumor","Metastatic" ))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
temp_data=Without_Treatment[which(Without_Treatment$Source == "Tumor"),]
data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]
#source_number=table(Without_Treatment$Source)/11
batch_calculate_signature_score<-function(temp_index){
  temp_data=Without_Treatment[which(Without_Treatment$Source == temp_index),]
  source_number=length(unique(temp_data$Sample))
  
  df3 <- data_summary(temp_data, varname="Proportion", 
                      groupnames=c("Cell_type"))
  df3$sd[which(df3$Proportion-df3$sd < 0)]=df3$Proportion[which(df3$Proportion-df3$sd < 0)]
  #df3$Cluster=as.factor(df3$Cluster)
  color_list=c("#DA5344","#85BC95","sandybrown","#9480B3","#2E78A5")
  color_index=match(temp_index,c("Blood","Normal","Precancerous","Tumor","Metastatic" ))
  df3$Cell_type=factor(df3$Cell_type,levels=Cell_order)
  p=ggplot(df3, aes(x=Cell_type, y=Proportion)) +
    geom_bar(stat="identity",color=color_list[color_index],fill=color_list[color_index], position=position_dodge(),alpha=0.9) +
    geom_errorbar(aes(ymin=Proportion-sd, ymax=Proportion+sd),color=color_list[color_index], width=.2,position=position_dodge(.9))+
    labs(title=paste0(temp_index," Sample number: ",source_number))+ theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 90))
  ggsave(file=paste0("Treat_naive_",temp_index,"_celltype_Proportion_bar.pdf"),p,width = 3.5,height = 2.5)
}
result=apply(matrix(levels(Without_Treatment$Source)),1,batch_calculate_signature_score)

#########For the treatment dataset, seperately calculated cell type proportion in pre and post treatment#####

Immthe_Treatment=sample_pro[which(sample_pro$Dataset %in% Immthe_Dataset),]
Immthe_Treatment$Treatment_timepoint=meta_info$Treatment_timepoint[match(Immthe_Treatment$Sample,meta_info$Sample)]
Immthe_Treatment=Immthe_Treatment[-1*which(Immthe_Treatment$Treatment_timepoint == "Chemo"),]
Immthe_Treatment$Source=factor(Immthe_Treatment$Source,levels=c("Blood","Normal","Tumor","Metastatic" ))
p=ggplot(Immthe_Treatment, aes(x=Cell_type, y=Proportion, color=Source,fill=Source)) + geom_boxplot(position=position_dodge(0.85),outlier.size =0,outlier.color ="white",alpha=0.9)+ theme_bw()+theme(axis.text.x = element_text(angle =90,vjust = 0.5),panel.background = element_rect( colour = "white"),legend.position = "top")+ scale_color_manual(values = c("#DA5344","#85BC95","#9480B3","#2E78A5"))+scale_fill_manual(values = c("#DA5344","#85BC95","#9480B3","#2E78A5"))+stat_compare_means( aes(label = ..p.signif..))
ggsave(file="TNK_distribution_across_condation_ImmTherapy.pdf",p,width = 12,height = 4)

temp_data=Immthe_Treatment[which(Immthe_Treatment$Source == "Tumor"),]
data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]

batch_calculate_signature_score<-function(temp_index){
  temp_data=Immthe_Treatment[which(Immthe_Treatment$Source == temp_index),]
  source_number=length(unique(temp_data$Sample))
  
  df3 <- data_summary(temp_data, varname="Proportion", 
                      groupnames=c("Cell_type"))
  df3$sd[which(df3$Proportion-df3$sd < 0)]=df3$Proportion[which(df3$Proportion-df3$sd < 0)]
  #df3$Cluster=as.factor(df3$Cluster)
  color_list=c("#DA5344","#85BC95","#9480B3","#2E78A5")
  color_index=match(temp_index,c("Blood","Normal","Tumor","Metastatic" ))
  df3$Cell_type=factor(df3$Cell_type,levels=Cell_order)
  p=ggplot(df3, aes(x=Cell_type, y=Proportion)) +
    geom_bar(stat="identity",color=color_list[color_index],fill=color_list[color_index], position=position_dodge(),alpha=0.9) +
    geom_errorbar(aes(ymin=Proportion-sd, ymax=Proportion+sd),color=color_list[color_index], width=.2,position=position_dodge(.9))+
    labs(title=paste0(temp_index," Sample number: ",source_number))+ theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 90))
  ggsave(file=paste0("Treat_Immther_",temp_index,"_celltype_Proportion_bar.pdf"),p,width = 3.5,height = 2.5)
}
result=apply(matrix(levels(Immthe_Treatment$Source)),1,batch_calculate_signature_score)


################Pre treatment###################

Immthe_Treatment$Treatment_timepoint=factor(Immthe_Treatment$Treatment_timepoint,levels=c("Pre-treatment","Post-treatment"))
p=ggplot(Immthe_Treatment, aes(x=Cell_type, y=Proportion, color=Treatment_timepoint,fill=Treatment_timepoint)) + geom_boxplot(position=position_dodge(0.85),outlier.size =0,outlier.color ="white",alpha=0.9)+ theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.background = element_rect( colour = "white"),legend.position = "top")+ scale_color_manual(values = c("#339B88","#4F6086"))+scale_fill_manual(values = c("#339B88","#4F6086"))+stat_compare_means( aes(label = ..p.signif..),method = "t.test")
ggsave(file="CD8TNK_distribution_across_condation_ImmTherapy_ttest.pdf",p,width = 6,height = 3.5)

temp_data=Immthe_Treatment[which(Immthe_Treatment$Treatment_timepoint == "Post-treatment"),]
data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]
batch_calculate_signature_score<-function(temp_index){
  temp_data=Immthe_Treatment[which(Immthe_Treatment$Treatment_timepoint == temp_index),]
  source_number=length(unique(temp_data$Sample))
  #data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
  #Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]
  df3 <- data_summary(temp_data, varname="Proportion", 
                      groupnames=c("Cell_type"))
  df3$sd[which(df3$Proportion-df3$sd < 0)]=df3$Proportion[which(df3$Proportion-df3$sd < 0)]
  #df3$Cluster=as.factor(df3$Cluster)
  color_list=c("#339B88","#4F6086")
  color_index=match(temp_index,c("Pre-treatment","Post-treatment"))
  df3$Cell_type=factor(df3$Cell_type,levels=Cell_order)
  p=ggplot(df3, aes(x=Cell_type, y=Proportion)) +
    geom_bar(stat="identity",color=color_list[color_index],fill=color_list[color_index], position=position_dodge(),alpha=0.9) +
    geom_errorbar(aes(ymin=Proportion-sd, ymax=Proportion+sd),color=color_list[color_index], width=.2,position=position_dodge(.9))+
    labs(title=paste0(temp_index," Sample number: ",source_number))+ theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 90))
  ggsave(file=paste0("Immthe_Treatment_",temp_index,"_celltype_Proportion_bar.pdf"),p,width = 3.5,height = 2.5)
}
result=apply(matrix(levels(Immthe_Treatment$Treatment_timepoint)),1,batch_calculate_signature_score)


#####Compare cell type proportion in Response and  non-response samples######
Respinse_meta=Immthe_Treatment[which(Immthe_Treatment$Treatment_timepoint == "Post-treatment"),]
Respinse_meta$Response=meta_info$Respponse[match(Respinse_meta$Sample,meta_info$Sample)]
Respinse_meta=Respinse_meta[-1*which(Respinse_meta$Response == "None"),]
Respinse_meta$Response[which(Respinse_meta$Response %in% c("MPR","R","Responded") )]="Responder"
Respinse_meta$Response[which(Respinse_meta$Response %in% c("Non-MPR","Non-responded","NR") )]="Non-responder"


Respinse_meta$Response=factor(Respinse_meta$Response,levels=c("Responder","Non-responder"))
p=ggplot(Respinse_meta, aes(x=Cell_type, y=Proportion, color=Response,fill=Response)) + geom_boxplot(position=position_dodge(0.85),outlier.size =0,outlier.color ="white",alpha=0.9)+ theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.background = element_rect( colour = "white"),legend.position = "top")+ scale_color_manual(values = c("#E64B35B2","#4DBBD5B2"))+scale_fill_manual(values = c("#E64B35B2","#4DBBD5B2"))+stat_compare_means( aes(label = ..p.signif..),method = "t.test")
ggsave(file="CD8TNK_distribution_across_condation_ImmTherapy_Respose_ttest.pdf",p,width = 6,height = 3.5)

temp_data=Respinse_meta[which(Respinse_meta$Response == "Non-responder"),]
data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]
batch_calculate_signature_score<-function(temp_index){
  temp_data=Respinse_meta[which(Respinse_meta$Response == temp_index),]
  source_number=length(unique(temp_data$Sample))
  #data_stat=aggregate(Proportion~Cell_type,temp_data,mean)
  #Cell_order=data_stat$Cell_type[order(data_stat$Proportion,decreasing = TRUE)]
  df3 <- data_summary(temp_data, varname="Proportion", 
                      groupnames=c("Cell_type"))
  df3$sd[which(df3$Proportion-df3$sd < 0)]=df3$Proportion[which(df3$Proportion-df3$sd < 0)]
  #df3$Cluster=as.factor(df3$Cluster)
  color_list=c("#E64B35B2","#4DBBD5B2")
  color_index=match(temp_index,c("Responder","Non-responder"))
  df3$Cell_type=factor(df3$Cell_type,levels=Cell_order)
  p=ggplot(df3, aes(x=Cell_type, y=Proportion)) +
    geom_bar(stat="identity",color=color_list[color_index],fill=color_list[color_index], position=position_dodge(),alpha=0.9) +
    geom_errorbar(aes(ymin=Proportion-sd, ymax=Proportion+sd),color=color_list[color_index], width=.2,position=position_dodge(.9))+
    labs(title=paste0(temp_index," Sample number: ",source_number))+ theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 90))
  ggsave(file=paste0("Immthe_Response_",temp_index,"_celltype_Proportion_bar.pdf"),p,width = 3.5,height = 2.5)
}
result=apply(matrix(levels(Respinse_meta$Response)),1,batch_calculate_signature_score)



############Proportition in cancer type#######
library(ggplot2)
library(RColorBrewer)
data=table(SeuratObj@meta.data$Cancer_type,SeuratObj@meta.data$curated_anno)
data=as.data.frame(data)
colnames(data)=c("CancerType","CellType","Fre")

batch_pie_for_each_celltype<-function(temp_celltype){
  temp_data=data[data$CellType == temp_celltype,]
  temp_data=droplevels(temp_data)
  temp_data$CancerType=paste(temp_data$CancerType,temp_data$Fre,sep="_")
  temp_data$CancerType=factor(temp_data$CancerType,levels=temp_data$CancerType[order(temp_data[,3],decreasing = TRUE)])
  
  
  # Barplot
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  bp<- ggplot(temp_data, aes(x="", y=Fre, fill=CancerType))+
    geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values =c(brewer.pal(n=12, name = "Paired"),brewer.pal(n=8, name = "Set2"),brewer.pal(n=12, name = "Set3"),brewer.pal(n=4, name = "Set1")))+blank_theme+theme(axis.text.x=element_blank())+ggtitle(temp_celltype)
  ggsave(file=paste0(temp_celltype,"_cancertype_dist.pdf"),bp,width = 6,height = 5)
  
}
result=apply(matrix(unique(data$CellType)),1,batch_pie_for_each_celltype)



stat_result=table(meta.data$Dataset,meta.data$curated_anno)

batch_cancer_type_propor<-function(temp_index){
  temp_result=stat_result[temp_index,]/sum(stat_result[temp_index,])
  return(temp_result)
}
result=apply(matrix(1:nrow(stat_result)),1,batch_cancer_type_propor)

colnames(result)=rownames(stat_result)

cutoff_1_result=result
cutoff_1_result[which(cutoff_1_result < 0.05)]=0
cutoff_1_result[which(cutoff_1_result >= 0.05)]=1
Convert_Dataset_Cancertype<-function(temp_cancer){
  mini_cancer=Cancer_type_class[which(Cancer_type_class$Major_Cancer == temp_cancer),1]
  Dataset=unique(meta.data$Dataset[which(meta.data$Cancer_type %in% mini_cancer)])
  if(length(Dataset) >1){
    return(rowSums(cutoff_1_result[,match(Dataset,colnames(cutoff_1_result))]))
  }else{
    return(cutoff_1_result[,match(Dataset,colnames(cutoff_1_result))])
  }
}
result=apply(matrix(unique(Cancer_type_class$Major_Cancer)),1,Convert_Dataset_Cancertype)

result=t(result)
rownames(result)=unique(Cancer_type_class$Major_Cancer)
#把健康的和cancer type的分开画，
#考虑到每个cancer type 中数据集的数目不一样多，而且比如 CD8_PDCD1_exh就真的只是在 NSCLC_GSE176021_aPD1 这套数据中出现了
library(pheatmap)
annotation_row=Cancer_type_class[match(rownames(result),Cancer_type_class[,2]),5]
annotation_row=as.data.frame(annotation_row)
rownames(annotation_row)=rownames(result)

result[which(result !=0 )]=1
library(ComplexHeatmap)
ha_row <- HeatmapAnnotation(Originated=annotation_row$annotation_row, which="row")

pdf("CancerType_Dataset5_Cell_proportion.pdf",width =6,height = 7 )
Heatmap(result,colorRampPalette(colors = c("white","#E64536"))(10), left_annotation =ha_row,row_split = data.frame(annotation_row))
dev.off()


############Propor in cancer type sepe ICB#####################
Datainfo=read.table("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/TISCH_V1V2_10X_data_patient_cell1.txt",header = TRUE,sep="\t")
Treatment_naive=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "None")]
Immthe_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Immunotherapy")]
Chemo_Dataset=Datainfo$Dataset.ID[ which(Datainfo$Treatment == "Chemotherapy")]

stat_result=table(meta.data$Dataset,meta.data$curated_anno)

batch_cancer_type_propor<-function(temp_index){
  temp_result=stat_result[temp_index,]/sum(stat_result[temp_index,])
  return(temp_result)
}
result=apply(matrix(1:nrow(stat_result)),1,batch_cancer_type_propor)

colnames(result)=rownames(stat_result)

cutoff_1_result=result
cutoff_1_result[which(cutoff_1_result < 0.01)]=0
cutoff_1_result[which(cutoff_1_result >= 0.01)]=1
Present_ICB_data<<-NULL
Convert_Dataset_Cancertype<-function(temp_cancer){
  mini_cancer=Cancer_type_class[which(Cancer_type_class$Major_Cancer == temp_cancer),1]
  Dataset=unique(meta.data$Dataset[which(meta.data$Cancer_type %in% mini_cancer)])
  ICB_data=intersect(Immthe_Dataset,Dataset)
  Non_ICB_data=setdiff(Dataset,ICB_data)
  
  if(length(ICB_data) >1){
    temp_icb=c(temp_cancer,"ICB_Treatment",length(which(meta.data$Dataset %in% ICB_data)),rowSums(cutoff_1_result[,match(ICB_data,colnames(cutoff_1_result))]))
  }
  if(length(ICB_data) == 1){
    temp_icb=c(temp_cancer,"ICB_Treatment",length(which(meta.data$Dataset %in% ICB_data)),cutoff_1_result[,match(ICB_data,colnames(cutoff_1_result))])
  }
  if(length(Non_ICB_data) >1){
    temp_nonicb=c(temp_cancer,"Treatment_naive",length(which(meta.data$Dataset %in% Non_ICB_data)),rowSums(cutoff_1_result[,match(Non_ICB_data,colnames(cutoff_1_result))]))
  }
  if(length(Non_ICB_data) == 1){
    temp_nonicb=c(temp_cancer,"Treatment_naive",length(which(meta.data$Dataset %in% Non_ICB_data)),cutoff_1_result[,match(Non_ICB_data,colnames(cutoff_1_result))])
  }
  
  if(length(Non_ICB_data) >= 1 && length(ICB_data)>=1){
    Present_ICB_data<<-rbind(Present_ICB_data,rbind(temp_icb,temp_nonicb))}
  if(length(Non_ICB_data) >= 1 && length(ICB_data)==0){
    Present_ICB_data<<-rbind(Present_ICB_data,temp_nonicb)}
  if(length(Non_ICB_data) == 0 && length(ICB_data)>=1 ){
    Present_ICB_data<<-rbind(Present_ICB_data,temp_icb)}
  
}
result=apply(matrix(unique(Cancer_type_class$Major_Cancer)),1,Convert_Dataset_Cancertype)

rownames(Present_ICB_data)=paste(Present_ICB_data[,1],Present_ICB_data[,2],sep="|")

#把健康的和cancer type的分开画，
#考虑到每个cancer type 中数据集的数目不一样多，而且比如 CD8_PDCD1_exh就真的只是在 NSCLC_GSE176021_aPD1 这套数据中出现了
library(pheatmap)
annotation_row=cbind(Cancer_type_class[match(Present_ICB_data[,1],Cancer_type_class[,2]),5],Present_ICB_data[,2])
annotation_row[which(annotation_row[,1] == "Healthy_derived" ),2]="Healthy_derived"
annotation_row=as.data.frame(annotation_row)
rownames(annotation_row)=rownames(Present_ICB_data)
#Cell_number_cancer_ICB=data.frame(Cell_number=as.numeric(Present_ICB_data[,3]))
Cell_number_cancer_ICB=data.frame(Cell_number=log2(as.numeric(Present_ICB_data[,3])))
rownames(Cell_number_cancer_ICB)=rownames(Present_ICB_data)

celltype_order=colnames(Present_ICB_data)[-1*c(1,2,3)]
Present_ICB_data=matrix(as.numeric(unlist(Present_ICB_data[,-1*c(1,2,3)])),byrow=FALSE,nrow=dim(Present_ICB_data)[1])
colnames(Present_ICB_data)=celltype_order
rownames(Present_ICB_data)=rownames(annotation_row)



library(ComplexHeatmap)

ha_row = rowAnnotation(ICB=annotation_row$V2,
                       Originated=annotation_row$V1,
                       Log_Cell_number = anno_barplot(Cell_number_cancer_ICB$Cell_number),
                       annotation_name_rot = 45
)


Present_ICB_data[which(Present_ICB_data !=0 )]=1
pdf("CancerType_Dataset_appear_Cell_proportion.pdf",width =8,height = 7 )
Heatmap(Present_ICB_data,colorRampPalette(colors = c("white","#E64536"))(10), left_annotation =ha_row,row_split = data.frame(annotation_row$V2))
dev.off()

#直接比在多少套数据中出现过，这个不公平。
#下面直接考虑用proportion， 每个cancer type在ICB的情况下

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

#把健康的和cancer type的分开画，
#考虑到每个cancer type 中数据集的数目不一样多，而且比如 CD8_PDCD1_exh就真的只是在 NSCLC_GSE176021_aPD1 这套数据中出现了
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
                       annotation_name_rot = 45
)


#Present_ICB_data[which(Present_ICB_data !=0 )]=1
pdf("CancerType_Dataset_Proportion_Cell_proportion.pdf",width =8,height = 7 )
Heatmap(Propor_ICB_data,c(colorRampPalette(colors = c("white","#E64536"))(20),rep("#E64536",30)), left_annotation =ha_row,row_split = data.frame(annotation_row$V2))
dev.off()