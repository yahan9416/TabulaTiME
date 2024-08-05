Mini_MPscore=readRDS("./Myeloid/Myeloid_Minicluster_Highqulaty_MPs_GSVA_score.rds")
Dataset_name=unlist(lapply(strsplit(colnames(Mini_MPscore),"@"),function(x) x[1]))

scMPscore_file=list.files("./ITH_NMF/Meta_program/Biology_vairation_SingleCellMetaCell/Myeloid/Batch_singlcell_level")
Metacell_info=list.files(".5/Mietacell/New_minicluster_exp_remDoublets")
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
  
  saveRDS(dataset_result,file=paste0("./ITH_NMF/Meta_program/Biology_vairation_SingleCellMetaCell/Myeloid/NewCorrelation_withinMoMacro/",dataset_name,"_singlecell_Metacell_cor.rds"))
}
result=apply(matrix(scMPscore_file[c(55:87)]),1,batch_process_scMetacell_MPscore_cor)
