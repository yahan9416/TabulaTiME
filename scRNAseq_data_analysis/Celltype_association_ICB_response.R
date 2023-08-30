#ICB response
#Collected the ICB response datasets from Tiger database

library(GSVA)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(MAESTRO)
library(Seurat)

Seurat_obj=
Seurat_obj@meta.data$curated_anno=droplevels(Seurat_obj@meta.data$curated_anno)
CD8TNK_DEG
CD8TNK_DEG=CD8TNK_DEG[order(CD8TNK_DEG$avg_log2FC,decreasing = TRUE),]
expmat=GetAssayData(Seurat_obj)

Tiger_ICB=list.files("/fs/home/hanya/Project/TME_Immune_difference/Tiger_ICB_Data")
ICB_dataset=gsub(".Rds","",Tiger_ICB)
ICB_dataset=gsub(".tsv","",ICB_dataset)
ICB_dataset=unique(ICB_dataset)


CD8TNK_DEG_cor_Mela<<-NULL
batch_investigate_ICB<-function(temp_icb_data){
  Melanoma=readRDS(paste0("./Tiger_ICB_Data/",temp_icb_data,".Rds"))
  Melanoma_cli=read.table(paste0("./Tiger_ICB_Data/",temp_icb_data,".tsv"),header = TRUE,sep="\t",fill = TRUE)
  gene_order=Melanoma[,1]
  Melanoma=as.matrix(Melanoma)
  Melanoma=Melanoma[,-1]
  Melanoma_exp=matrix(as.numeric(unlist(Melanoma)),byrow=FALSE,ncol=dim(Melanoma)[2])
  colnames(Melanoma_exp)=colnames(Melanoma)
  rownames(Melanoma_exp)=as.vector(unlist(gene_order))
  Melanoma_exp[which(is.na(Melanoma_exp))]=0
  #Remove genes which are expression level are 0 in all samples 
  Melanoma_exp=Melanoma_exp[-1*which(rowSums(Melanoma_exp) == 0),]
  
  
  batch_subtract_mean<-function(x){
    return(c(Melanoma_exp[x,]-mean(Melanoma_exp[x,])))
  }
  coore_Mealnoma_exp=apply(matrix(1:nrow(Melanoma_exp)),1,batch_subtract_mean)
  coore_Mealnoma_exp=t(coore_Mealnoma_exp)
  rownames(coore_Mealnoma_exp)=rownames(Melanoma_exp)
  if(temp_icb_data %in% ICB_dataset[6]){
    Melanoma_cli$Treatment[which( Melanoma_cli$Treatment %in%  c("POST1","POST2","POST3","POST4","POST5"))]="POST"
  }
  treatment_response_sta=table(Melanoma_cli$Treatment,Melanoma_cli$response_NR)
  
  batch_select_celltype_DEG<-function(temp_celltype){
    temp_DEG_gene=CD8TNK_DEG[which(CD8TNK_DEG$cluster == temp_celltype),]
    DEG_gene=na.omit(intersect(temp_DEG_gene$gene,rownames(coore_Mealnoma_exp))[1:100])
    return(DEG_gene)
  }
  DEG_gene_list=apply(matrix(unique(CD8TNK_DEG$cluster)),1,batch_select_celltype_DEG)
  names(DEG_gene_list)=unique(CD8TNK_DEG$cluster)
  
  #EDT early during treatment
  #return(table(Melanoma_cli$Treatment))
  #这里面可以分为pre, on ,post 三种，要分开吗？分开吧，
  #没有的我们暂且把它当做pre，然后我们要求每个分组里面至少有两个sample
  if(length(table(Melanoma_cli$Treatment)) == 0){
    DEG_GSVA=gsva(coore_Mealnoma_exp,DEG_gene_list)
    
    batch_cell_type_ICB_response<-function(temp_celltype){
      
      temp_data=cbind(DEG_GSVA[match(temp_celltype,rownames(DEG_GSVA)),],colnames(coore_Mealnoma_exp),temp_celltype,Melanoma_cli$response_NR[match(colnames(coore_Mealnoma_exp),Melanoma_cli$sample_id)])
      temp_data=as.data.frame(temp_data)
      colnames(temp_data)=c("GSVA","Sample","Celltype","Response")
      temp_data$GSVA=as.numeric(as.vector(temp_data$GSVA))
      temp_stat=c(temp_celltype,temp_icb_data,"WithoutGroup",unlist(table(temp_data$Response)),unlist(aggregate(GSVA~Response,temp_data,mean)),unlist(wilcox.test(GSVA~Response,temp_data))[2])
      CD8TNK_DEG_cor_Mela<<-rbind(CD8TNK_DEG_cor_Mela,temp_stat)
    } 
    result=apply(matrix(unique(CD8TNK_DEG$cluster)),1,batch_cell_type_ICB_response)
    
  }
  
  if(length(table(Melanoma_cli$Treatment)) > 0){
    rownames(treatment_response_sta)[which(rownames(treatment_response_sta) == "")]="WithoutGroup"
    Melanoma_cli$Treatment[which(Melanoma_cli$Treatment == "")]="WithoutGroup"
    batch_subgroup<-function(temp_group){
      if(!is.na(match("UNK",colnames(treatment_response_sta)))){
        temp_treatment_response_sta=treatment_response_sta[temp_group,-1*match("UNK",colnames(treatment_response_sta))]
        min_group_sample=min(temp_treatment_response_sta)
      }else{
        min_group_sample=min(treatment_response_sta[temp_group,])
      }
      if(min_group_sample >= 2){
        
        coore_Mealnoma_exp_group=coore_Mealnoma_exp[,na.omit(match(Melanoma_cli$sample_id[which(Melanoma_cli$Treatment == temp_group)],colnames(coore_Mealnoma_exp)))]
        DEG_GSVA=gsva(coore_Mealnoma_exp_group,DEG_gene_list)
        
        batch_subgroup_cell_type_ICB_response<-function(temp_celltype){
          
          
          temp_data=cbind(DEG_GSVA[match(temp_celltype,rownames(DEG_GSVA)),],colnames(coore_Mealnoma_exp_group),temp_celltype,Melanoma_cli$response_NR[match(colnames(coore_Mealnoma_exp_group),Melanoma_cli$sample_id)])
          temp_data=as.data.frame(temp_data)
          colnames(temp_data)=c("GSVA","Sample","Celltype","Response")
          if("UNK" %in% temp_data$Response){
            temp_data=temp_data[-1*(which(temp_data$Response == "UNK" )),]
          }
          temp_data$GSVA=as.numeric(as.vector(temp_data$GSVA))
          temp_stat=c(temp_celltype,temp_icb_data,temp_group,unlist(table(temp_data$Response)),unlist(aggregate(GSVA~Response,temp_data,mean)),unlist(wilcox.test(GSVA~Response,temp_data))[2])
          CD8TNK_DEG_cor_Mela<<-rbind(CD8TNK_DEG_cor_Mela,temp_stat)
        } 
        result=apply(matrix(unique(CD8TNK_DEG$cluster)),1,batch_subgroup_cell_type_ICB_response)
      }
      
    }
    temp_result=apply(matrix(rownames(treatment_response_sta)),1,batch_subgroup)
  } 
}
result=apply(matrix(ICB_dataset),1,batch_investigate_ICB)

#result=apply(matrix(ICB_dataset[15:27]),1,batch_investigate_ICB)



colnames(CD8TNK_DEG_cor_Mela)[1:5]=c("Celltype","Dataset","Group","N_numnber","R_number")
CD8TNK_DEG_cor_Mela=as.data.frame(CD8TNK_DEG_cor_Mela)
CD8TNK_DEG_cor_Mela$Correlation1=as.numeric(CD8TNK_DEG_cor_Mela$Correlation1)
CD8TNK_DEG_cor_Mela$Correlation2=as.numeric(CD8TNK_DEG_cor_Mela$Correlation2)
CD8TNK_DEG_cor_Mela$p.value=as.numeric(CD8TNK_DEG_cor_Mela$p.value)
CD8TNK_DEG_cor_Mela$Direction=CD8TNK_DEG_cor_Mela$Correlation2 > CD8TNK_DEG_cor_Mela$Correlation1

batch_get_gap_between_RNR<-function(temp_index){
  temp_data=CD8TNK_DEG_cor_Mela[temp_index,]
  if(length(which(temp_data[8:9] < 0 )) == 1){
    gap=sum(abs(temp_data[8:9]))}
  if(length(which(temp_data[8:9] < 0 )) == 2){
    gap=abs(abs(temp_data[8])-abs(temp_data[9]))}
  if(length(which(temp_data[8:9] < 0 )) == 0){
    gap=abs(temp_data[8]-temp_data[9])}
  return(gap)
}
result=apply(matrix(1:dim(CD8TNK_DEG_cor_Mela)[1]),1,batch_get_gap_between_RNR)
CD8TNK_DEG_cor_Mela$gap=as.vector(unlist(result))
CD8TNK_DEG_cor_Mela$Direction=ifelse(CD8TNK_DEG_cor_Mela$Direction, 1, -1)
CD8TNK_DEG_cor_Mela$gap=CD8TNK_DEG_cor_Mela$Direction * CD8TNK_DEG_cor_Mela$gap


gap=matrix(as.vector(CD8TNK_DEG_cor_Mela$gap),byrow=TRUE,ncol = length(unique(CD8TNK_DEG_cor_Mela[,1])) )
colnames(gap)=unique(CD8TNK_DEG_cor_Mela$Celltype)
rownames(gap)=unique(paste(CD8TNK_DEG_cor_Mela$Dataset,CD8TNK_DEG_cor_Mela$Group,sep="|"))
p_value=matrix(as.vector(CD8TNK_DEG_cor_Mela$p.value),byrow=TRUE,ncol = length(unique(CD8TNK_DEG_cor_Mela[,1]))  )
colnames(p_value)=unique(CD8TNK_DEG_cor_Mela$Celltype)
rownames(p_value)=unique(paste(CD8TNK_DEG_cor_Mela$Dataset,CD8TNK_DEG_cor_Mela$Group,sep="|"))

if (!is.null(p_value)){
  ssmt <- p_value< 0.01
  p_value[ssmt] <-'**'
  smt <- p_value >0.01& p_value <0.05
  p_value[smt] <- '*'
  p_value[!ssmt&!smt]<- ''
} else {
  p_value <- F
}
annotation=data.frame(group=unlist(lapply(strsplit(rownames(gap),"\\|"),function(x) x[2])))
annotation$group[which(annotation$group == "EDT")]="ON"
annotation$group[which(annotation$group == "Normal")]="PRE"



library(ComplexHeatmap)
ha_row = rowAnnotation(Treatment=annotation$group,
                       col = list(Treatment=c("WithoutGroup"="#91e7fc","PRE"="#ccb6e4","POST"="#b5fe8e","ON"="#aed2e1")))
pdf("Fibroblast_Top100_corre_Tiger_Spe_add_PrePost_complex_add_removeGBM.pdf",width =7,height = 7 )
bk <- c(seq(-1,-0.1,by=0.1),seq(0.1,1,by=0.1))
Heatmap(gap,c(rep("#E64536",8),colorRampPalette(colors = c("#E64536","white"))(length(bk)/2),colorRampPalette(colors = c("white","#414393"))(length(bk)/2),rep("#414393",8)), row_split = data.frame(annotation$group),left_annotation =ha_row,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text( p_value[i, j], x, y, gp = gpar(fontsize = 15, col = "white"))
        })
dev.off()

