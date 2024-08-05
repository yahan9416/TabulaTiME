library(mclust)
#1-ARIbatch


set.seed(0)
cluster_label=Seurat_obj@meta.data$seurat_clusters
cluster_batch=Seurat_obj@meta.data$batch
batch_random_20<-function(temp_index){
  temp_index=sample(1:length(cluster_batch),5000)
  ari_batch <- adjustedRandIndex(cluster_batch[temp_index], cluster_label[temp_index])
  return(ari_batch)
}
result=apply(matrix(1:20),1,batch_random_20)

#ARIcelltype
cluster_label=Seurat_obj@meta.data$seurat_clusters
cluster_batch=Seurat_obj@meta.data$annotation
batch_random_20<-function(temp_index){
  temp_index=sample(1:length(cluster_batch),5000)
  ari_celltype <- adjustedRandIndex(cluster_batch[temp_index], cluster_label[temp_index])
  return(ari_celltype)
}
result=apply(matrix(1:20),1,batch_random_20)
