# Merge all minicluster(metacell)


Minicluster_file=list.files("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Minicluster")
Minicluster_file=Minicluster_file[grep("RData",Minicluster_file)]
#First: Count how many datasets each gene appeared in
setwd("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Minicluster")
Minicluster_TPM_gene<<-NULL
batch_merge_mini_cluster<-function(temp_index){
  load(temp_index)
  rownames(Minicluster_TPM)=gsub('\\"',"",rownames(Minicluster_TPM))
  return(rownames(Minicluster_TPM))
}
Minicluster_TPM_gene=apply(matrix(Minicluster_file),1,batch_merge_mini_cluster)
length(which(table(as.vector(unlist(Minicluster_TPM_gene))) >= 60 ))

#And then selected the gene used to 
commone_gene=names(table(as.vector(unlist(Minicluster_TPM_gene))))[which(table(as.vector(unlist(Minicluster_TPM_gene))) >= 60 )]

Minicluster_TPM_gene=unique(as.vector(unlist(Minicluster_TPM_gene))) 

Minicluster_file=list.files("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Minicluster")
Minicluster_TPM_Total_1<<-NULL
batch_merge_mini_cluster<-function(temp_index){
  print(temp_index)
  load(temp_index)
  rownames(Minicluster_TPM)=gsub('\\"',"",rownames(Minicluster_TPM))
  #index=na.omit(match(rownames(Minicluster_TPM),Minicluster_TPM_gene))
  inter_gene=intersect(rownames(Minicluster_TPM),Minicluster_TPM_gene)
  total_gene_index=match(inter_gene,Minicluster_TPM_gene)
  temp_gene_index=match(inter_gene,rownames(Minicluster_TPM))
  
  temp_tpm=matrix(0, ncol=dim(Minicluster_TPM)[2],nrow = length(Minicluster_TPM_gene) )
  temp_tpm[total_gene_index,]=Minicluster_TPM[temp_gene_index,]
  
  #temp_data=unlist(strsplit(minicluster_file[temp_index],"/"))[11]
  temp_data=gsub("_Minicluster_Data.RData","",temp_index)
  colnames(temp_tpm)=paste(temp_data,colnames(Minicluster_TPM),sep="@")
  Minicluster_TPM_Total_1<<-cbind(Minicluster_TPM_Total_1,temp_tpm)
  #return(temp_tpm)
}
result=apply(matrix(Minicluster_file),1,batch_merge_mini_cluster)
rownames(Minicluster_TPM_Total_1)<-Minicluster_TPM_gene

