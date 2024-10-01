#' Calculating the correlation of the MP signature score between metacells and their corresponding cells
#' @param Metacell_MPscore_path Path of calculated MP score at metacell level
#' @param singcell_MPscore_path Path of calculated MP score at single-cell level
#' @param Metacell_info_path Path of metacell information
#' @author Ya Han
#' 
args <- commandArgs(trailingOnly = TRUE)
Metacell_MPscore_path <- args[1] #The path of MP signature score of each Metacell
singcell_MPscore_path <- args[2]
Metacell_info_path <- args[3]


Mini_MPscore=readRDS(Metacell_MPscore_path)
Dataset_name=unlist(lapply(strsplit(colnames(Mini_MPscore),"@"),function(x) x[1]))

scMPscore_file=list.files(singcell_MPscore_path)
Metacell_info=list.files(Metacell_info_path)
MetaSc_MP_cor<<-NULL
batch_process_scMetacell_MPscore_cor<-function(temp_dataset){
  print(temp_dataset)
  scSignature_score=readRDS(temp_dataset)
  dataset_name=gsub("_sclevel_Highqulaty_MPs_GSVA_score.rds","",temp_dataset)
  load(Metacell_info[grep(dataset_name,Metacell_info)])
  index=which(Dataset_name == dataset_name) 
  colname_source_cluster=colnames(Minicluster_TPM)
  colname_source_cluster=paste0(dataset_name,"@",colname_source_cluster)
  Fibroblast_metacell_list=rownames(temp_meta)[which(temp_meta$assign.curated == "Mono/Macro")]
  
  batch_calculcate_MPscore_correlation_MiniSc<-function(temp_index){
    #print(temp_index)
    aa=unlist(strsplit(colnames(Mini_MPscore)[temp_index],"@"))[2]
    if(aa %in% Fibroblast_metacell_list){  
      #print(temp_index)
      Mini_index=match(colnames(Mini_MPscore)[temp_index],colname_source_cluster)
      if(exists("Cell_minicluster_list_new")){
        Mini_cells=Cell_minicluster_list_new[[Mini_index]]
      }else{
        Mini_cells=Cell_minicluster_list[[Mini_index]]}
      print(length(na.omit(match(Mini_cells,colnames(scSignature_score)))))
      if(length(na.omit(match(Mini_cells,colnames(scSignature_score)))) >0){
        temp_MPscore=scSignature_score[,match(Mini_cells,colnames(scSignature_score))]
        cor_reuslt=apply(temp_MPscore,2,function(x){return(unlist(cor.test(Mini_MPscore[,temp_index],x))[4])})
        return(cor_reuslt)}
    }
    
  }
  dataset_result=apply(matrix(index),1,batch_calculcate_MPscore_correlation_MiniSc)
  
  saveRDS(dataset_result,file=paste0(dataset_name,"_singlecell_Metacell_cor.rds"))
}
result=apply(matrix(scMPscore_file),1,batch_process_scMetacell_MPscore_cor)
