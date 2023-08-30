#Look at the immune checkpoint response effect of cell type upregulated gene on immune therapy dataset
#Using the ratio of specific cell type compare to all cells expression value to look at the effect on immune checkpoint therapy
#' @param Seurat_obj    The input is a seurat obj, the meta information of seurat need including curated_anno.
#' @param upDEG_list    The cell-type/cluster up-regulated gene list
#' @param ICB_exp       The path of expression of patient which under immune therapy.
#' @param Response_info The file of clinical information including the response or non-response
#' @param output_dir    The path of plot file store


library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(org.Hs.eg.db)

Correlation_ICB<-function(Seurat_obj,upDEG_list=NULL,ICB_exp,Response_info,output_dir=NULL){
  if(length(upDEG_list) == 0){
    upDEG_list=FindAllMarkers(Seurat_obj,only.pos = TRUE)
  }
  if(length(output_dir) == 0){output_dir=getwd()}
  
  
  #用在all TNK 的cell中找打的DEG,但是backgroud 是所有的Treg 细胞，
  expmat=GetAssayData(Seurat_obj)
  expmat_gene=rownames(expmat)
  expmat_rowname <- bitr(expmat_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  expmat=as.matrix(expmat)
  
  data_icb_ratio_all<<-NULL
  patient_order=NULL
  patient_res=NULL
  ICB_treatfiles=list.files(ICB_exp)  
  batch_calculate_response_ICB_file<-function(x){
    temp_list=list.files(paste0(ICB_exp,x))
    ICB_expMat=read.table(paste0(ICB_exp,x,"/",temp_list[2]),header = T,sep="\t",row.names = 1)
    Response=read.table(paste0(ICB_exp,x,"/",temp_list[1]),header = T,sep="\t")
    colnames(Response)=gsub("X","patient",colnames(Response))
    colnames(Response)=gsub("Patient","patient",colnames(Response))
    if(length(Response$patient) == 0){
      Response=cbind(rownames(Response),Response)
      colnames(Response)[1]="patient"
    }
    
    data_icb=NULL
    get_correlation_ratio<-function(clus){
      single_gene <- bitr(cluster_DEG$gene[which(cluster_DEG$cluster == clus)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      gene_intersect=intersect(single_gene$ENTREZID,rownames(ICB_expMat))
      
      if(length(Response$patient) >0  ){
        if(length(grep("-",Response$patient))>0){
          Response$patient=gsub("-","\\.",Response$patient)}
        colnames(ICB_expMat)=gsub("^X","",colnames(ICB_expMat))
        patient_res=Response$Response[match(colnames(ICB_expMat),Response$patient)]
        patient_order=Response$patient[match(colnames(ICB_expMat),Response$patient)]
        patient_res[which(patient_res == 0)]="NR"
        patient_res[which(patient_res == 1)]="R"
        ICB_expMat=as.matrix(ICB_expMat)#must be a matrix
        
        index_icb=match(gene_intersect,rownames(ICB_expMat))
        index_seu=match(expmat_rowname$SYMBOL[match(gene_intersect,expmat_rowname$ENTREZID)],expmat_gene)
        Treg_average=rowMeans(as.matrix(expmat))
        cluster_average=rowMeans(as.matrix(expmat[,which(TNK_Seurat@meta.data$curated_anno_gene == clus)]))
        cluster_average_ratio=(cluster_average/Treg_average)[index_seu]
        cor_result=apply(ICB_expMat[index_icb,],2,function(x){unlist(cor.test(x,cluster_average_ratio))[4]})
        temp_data=as.numeric(cor_result)
        data_icb<<-rbind(data_icb,cbind(temp_data,patient_order,patient_res,clus))
      }
    }
    result=apply(matrix(unique(cluster_DEG$cluster)),1,get_correlation_ratio)
    
    colnames(data_icb)[1]="Correlation"
    data_icb=as.data.frame(data_icb)
    data_icb$Correlation=as.numeric(as.vector(unlist(data_icb$Correlation)))
    
    setwd(output_dir)
    my_comparisons <- list( c("NR", "R"))
    p=ggplot(data_icb, aes(x=clus, y=Correlation,fill=patient_res)) + 
      geom_boxplot()+
      theme_classic()+ scale_fill_brewer(palette="Blues") + theme_classic()+              theme(legend.position="bottom")+theme(axis.text.x = element_text(angle = 45))+ stat_compare_means(label = "p.signif", paired = TRUE)
    ggsave(paste0("icb_response_",x,"_spear_.png"),p)
    
  }
  result=apply(as.matrix(ICB_treatfiles),1,batch_calculate_response_ICB_file)
  
}
