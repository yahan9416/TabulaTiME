minicluster_file=list.files("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Test_miniCluster_Performance/Minicluster",full.names = TRUE,recursive = TRUE)
batch_check<-function(temp_data){
  temp_minicluster_file=minicluster_file[grep(paste0(temp_data,"_Minicluster_"),minicluster_file)]
  #return(length(temp_minicluster_file))
  Coverage_mini<<-NULL
  batch_gene_coverage_eblow<-function(temp_file){
    load(temp_file)
    test_gene_coverage<-function(sample_gene_exp){
      return(length(which(sample_gene_exp != 0)))
    }
    gene_coverage=apply(Minicluster_TPM,2,test_gene_coverage)
    Coverage=data.frame(coverage=gene_coverage,cell_number=unlist(strsplit(temp_file,"\\/"))[10])
    Coverage_mini<<-rbind(Coverage_mini,Coverage)
  }
  result=apply(matrix(temp_minicluster_file),1,batch_gene_coverage_eblow)
  
  p=ggplot(Coverage_mini, aes(x=cell_number, y=coverage))+ggtitle(temp_data)+
    geom_boxplot( fill="#276AA7",color="black",alpha=0.8,outlier.size=0)+
    theme_classic()+theme(axis.text.x = element_text(angle = 90))+scale_x_discrete(limits=unique(Coverage_mini$cell_number)[c(1,5:12,2,3,4)])
  ggsave(paste0(temp_data,"_gene_coverage.pdf"),p,width=5,height=4)
  
  #################Mean#################
  CHOL_GINI_mean=aggregate(coverage~cell_number,Coverage_mini,mean)
  p=ggplot(data=CHOL_GINI_mean, aes(x=cell_number, y=coverage, group=1)) +
    geom_line()+theme(axis.text.x = element_text(angle = 90))+ggtitle(temp_data)+
    geom_point()+scale_x_discrete(limits=unique(CHOL_GINI_mean$cell_number)[c(1,5:12,2,3,4)])
  ggsave(paste0(temp_data,"_gene_coverage_mean.pdf"),p)
  
  rownames(CHOL_GINI_mean)=1:12
  CHOL_GINI_mean$cell_number=gsub("Minicluster","",CHOL_GINI_mean$cell_number)
  CHOL_GINI_mean$cell_number=lapply(strsplit(CHOL_GINI_mean$cell_number,"_"),function(x) x[1])
  CHOL_GINI_mean$cell_number=as.numeric(as.vector(CHOL_GINI_mean$cell_number))
  CHOL_GINI_mean=CHOL_GINI_mean[order(CHOL_GINI_mean$cell_number),]
  CHOL_GINI_mean$coverage=CHOL_GINI_mean$coverage
  
  aa=elbow_finder(CHOL_GINI_mean$cell_number,CHOL_GINI_mean$coverage)
  double_derived_point=which(aa$d2 == max(aa$d2,na.rm = TRUE))
  
  CHOL_GINI_mean=CHOL_GINI_mean[1:11,]
  aa=elbow_finder(CHOL_GINI_mean$cell_number,CHOL_GINI_mean$coverage)
  double_derived_point1=which(aa$d2 == max(aa$d2,na.rm = TRUE))
  #elbow_point=elbow.point(CHOL_GINI_mean$coverage)
  
  #################Median#################
  #Caluclate eblow for median of gene coverage
  CHOL_GINI_mean=aggregate(coverage~cell_number,Coverage_mini,median)
  p=ggplot(data=CHOL_GINI_mean, aes(x=cell_number, y=coverage, group=1)) +
    geom_line()+theme(axis.text.x = element_text(angle = 90))+ggtitle(temp_data)+
    geom_point()+scale_x_discrete(limits=unique(CHOL_GINI_mean$cell_number)[c(1,5:12,2,3,4)])
  ggsave(paste0(temp_data,"_gene_coverage_medain.pdf"),p)
  
  rownames(CHOL_GINI_mean)=1:12
  CHOL_GINI_mean$cell_number=gsub("Minicluster","",CHOL_GINI_mean$cell_number)
  CHOL_GINI_mean$cell_number=lapply(strsplit(CHOL_GINI_mean$cell_number,"_"),function(x) x[1])
  CHOL_GINI_mean$cell_number=as.numeric(as.vector(CHOL_GINI_mean$cell_number))
  CHOL_GINI_mean=CHOL_GINI_mean[order(CHOL_GINI_mean$cell_number),]
  
  aa=elbow_finder(CHOL_GINI_mean$cell_number,CHOL_GINI_mean$coverage)
  double_derived_point_med=which(aa$d2 == max(aa$d2,na.rm = TRUE))
  
  CHOL_GINI_mean=CHOL_GINI_mean[1:11,]
  aa=elbow_finder(CHOL_GINI_mean$cell_number,CHOL_GINI_mean$coverage)
  double_derived_point1_med=which(aa$d2 == max(aa$d2,na.rm = TRUE))
  
  return(c(temp_data,double_derived_point,double_derived_point1,double_derived_point_med,double_derived_point1_med))
  
}
result=apply(matrix(Select_dataset),1,batch_check)
result=t(result)
result=as.data.frame(result)
colnames(result)=c("Dataset","Mean_Eblow_point_all","Mean_Eblow_point_1_110","Median_Eblow_point_all","Median_Eblow_point_1_110")
