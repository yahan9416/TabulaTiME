# function
library(MAESTRO)
library(Seurat)
library(ggplot2)

#这个函数没有加上cell type和source的信息，希望可以一次性获得
Generate_DiffCell_minicluster_for_Seurat<-function(Seurat_obj,CellNumber){
  
  TNK_Seurat=Seurat_obj
  cluster_sample_average<-function(cluster){
    
    sample_index=rownames(TNK_Seurat@meta.data)[which( TNK_Seurat@meta.data$sample_cluster ==cluster)]
    
    if(length(sample_index) > CellNumber){
      temp_seurat=subset(TNK_Seurat,cells=sample_index)
      
      cells <- length(sample_index)
      cluster.res=0.2
      dims.use=1:15
      temp_seurat <- FindNeighbors(temp_seurat, dims = dims.use)
      temp_seurat=FindClusters(temp_seurat,resolution = cluster.res)
      #print(cluster)
      #print(table(temp_seurat@meta.data$seurat_clusters))
      temp_seurat@meta.data$source_cluster=temp_seurat@meta.data$seurat_clusters
      
      Find_the_near_100Neig_incluster<-function(y){
        print(y)
        umap_df <<- temp_seurat@reductions$umap@cell.embeddings[which(temp_seurat@meta.data$source_cluster == y),]
        cell_dist<<-as.matrix(dist(umap_df))
        if(nrow(cell_dist) > CellNumber){
          #if more than 50 cell will generate a new cluster or will not
          if(ceiling(length(which(temp_seurat@meta.data$source_cluster == y))/CellNumber) == round(length(which(temp_seurat@meta.data$source_cluster == y))/CellNumber) ){
            temp_cluster<<-ceiling(length(which(temp_seurat@meta.data$source_cluster == y))/CellNumber)
            #first generate the first cluster and remove less than 100 cells from the total, because there less than temp_cluster*100, so there must with overlap between clusters
            seed_cell=sample(rownames(cell_dist),1)
            cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:CellNumber])
            Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
            
            colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sep="|"))
            less10_cell_number=(length(which(temp_seurat@meta.data$source_cluster == y)) - floor(length(which(temp_seurat@meta.data$source_cluster == y))/CellNumber)*CellNumber)
            
            
            index=match(cluster_cells,rownames(cell_dist))
            cell_dist<<-cell_dist[-1*index,-1*index]
            
            generate_50cell_expresion<-function(sub_cluster){
              
              if(sub_cluster != temp_cluster){
                seed_cell=sample(rownames(cell_dist),1)
                cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:CellNumber])
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                index=match(cluster_cells,rownames(cell_dist))
                cell_dist<<-cell_dist[-1*index,-1*index]
              }else{
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(rownames(cell_dist)))
              }
              colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
            } 
            result=apply(as.matrix(2:temp_cluster),1,generate_50cell_expresion)
            
          }else{
            temp_cluster<<-round(length(which(temp_seurat@meta.data$source_cluster == y))/CellNumber)
            generate_50cell_expresion<-function(sub_cluster){
              seed_cell=sample(rownames(cell_dist),1)
              if(sub_cluster==temp_cluster){
                cluster_cells=colnames(cell_dist)[1:CellNumber]
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
              }else{
                cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:CellNumber])
                cell_index=match(cluster_cells,colnames(cell_dist))
                cell_dist<<-cell_dist[-1*cell_index,-1*cell_index]
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
              }
              
            } 
            result=apply(as.matrix(1:temp_cluster),1,generate_50cell_expresion)
          }
          
        }else{
          #cluster including total less than 100 cells 
          
          Cell_minicluster_list<<-c(Cell_minicluster_list,list(colnames(cell_dist)))
          colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sep="|"))
        }
        
      } 
      #source_cluster_moerthan100=names(which(table(temp_seurat@meta.data$source_cluster) >= 100))
      result=apply(as.matrix(unique(temp_seurat@meta.data$seurat_clusters)),1,Find_the_near_100Neig_incluster)
      
    }else{
      
      Cell_minicluster_list<<-c(Cell_minicluster_list,list(sample_index))
      
      colname_source_cluster<<-c(colname_source_cluster,paste(cluster,0,sep="|"))
      
    }
    
  }
  result=apply(as.matrix(unique(as.vector(unlist( TNK_Seurat@meta.data$sample_cluster)))),1,cluster_sample_average)   
  
}

#application
TPM=readRDS("./NSCLC_GSE176021_TPM_exp.rds")
Seurat_obj=readRDS("./NSCLC_GSE176021_aPD1_res.rds")
rownames(Seurat_obj@meta.data)=names(Idents(Seurat_obj))
if(length(Seurat_obj@meta.data$Sample) == 0){
  Seurat_obj@meta.data$sample_cluster=paste(Seurat_obj@meta.data$Patient,Seurat_obj@meta.data$seurat_clusters,sep = "|")
}else{
  Seurat_obj@meta.data$sample_cluster=paste(Seurat_obj@meta.data$Sample,Seurat_obj@meta.data$seurat_clusters,sep = "|")}


#cluster_sample_average 这个函数是用来产mini-cluster的，但是需要注意的是，我们现在要用不同的细胞数目
Cell_minicluster_list<<-list()
colname_source_cluster<<-NULL
Generate_DiffCell_minicluster_for_Seurat(Seurat_obj,30)

