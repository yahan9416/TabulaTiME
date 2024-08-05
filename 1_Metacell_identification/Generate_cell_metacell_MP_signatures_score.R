library(MAESTRO)
library(Seurat)
library(GSVA)

################Single cell######################
Myeloid_seu=readRDS("./Merge_minicluster30_Myeloid_Seurat_CCA_kidney_30.rds")
#including all MPs gene list after quality control
MPs_genelist=readRDS("./Myeloid_PotFoun_cluster_MPs_genelist.rds")#including all MPs gene list without quality control


batch_calculate_MPs_score_singlecell_level<-function(temp_dataset){
  temp_filelist=list.files(paste0("./TISCH1/static/data/",temp_dataset[2]),full.names = TRUE)
  Myeloid_seu=readRDS(temp_filelist[grep("_res.rds",temp_filelist)])
  Myeloid_seu=Myeloid_seu$RNA
  myofibro_cells=rownames(Myeloid_seu@meta.data)[which(Myeloid_seu@meta.data$assign.curated %in% c("Mono/Macro"))]
  Myeloid_seu=subset(Myeloid_seu,cells=myofibro_cells)
  DefaultAssay(Myeloid_seu)<-"RNA"
  expmat=GetAssayData(Myeloid_seu)
  expmat=as.matrix(expmat)
  
  Myeloid_Minicluster_MPs_score=gsva(expmat,MPs_genelist,method="ssgsea")
  saveRDS(Myeloid_Minicluster_MPs_score,file=paste0("./Myeloid/Batch_singlcell_level/",temp_dataset[1],"_sclevel_Highqulaty_MPs_GSVA_score.rds"))
}
result=apply(as.matrix(Myeloid_invo_dataset_tumor[5,]),1,batch_calculate_MPs_score_singlecell_level)



###################metacell level####################

Myeloid_seu=readRDS("./PanCancer_Minicluster30_Myelid_SeuratObj.rds")
Myeloid_seu@meta.data$curated_anno=as.vector(unlist(Myeloid_meta$curated_anno[match(rownames(Myeloid_seu@meta.data),rownames(Myeloid_meta))]))
#including all MPs gene list after quality control
MPs_genelist=readRDS("./ITH_NMF/Meta_program/Myeloid/Myeloid_PotFoun_cluster_MPs_genelist.rds")

DefaultAssay(Myeloid_seu)<-"RNA"
expmat=GetAssayData(Myeloid_seu)
expmat=as.matrix(expmat)

Myeloid_Minicluster_MPs_score=gsva(expmat,MPs_genelist,method="ssgsea")


