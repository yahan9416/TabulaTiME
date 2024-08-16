library(MAESTRO)
library(Seurat)
library(Biobase)
library(NMF)

#Generate intratumor heterogeneity by NMF for each sample
file_list=list.files("./data")
batch_process_nmf_program<-function(temp_dataset){
  temp_filelist=list.files(paste0("/fs/home/wangyuting/biosoft/django/wyt_test/TISCH1/static/data/",temp_dataset),full.names = TRUE)
  Seuratobj=readRDS(temp_filelist[grep("_res.rds",temp_filelist)])
  Seuratobj=Seuratobj$RNA
 
  Seuratobj@meta.data$assign.curated[which(Seuratobj@meta.data$assign.curated %in% c("CD4Tconv","CD8T","CD8Tex","Tprolif","Treg"))]="Tcell"
  Seuratobj@meta.data$assign.curated[which(Seuratobj@meta.data$assign.curated %in% c("DC","pDC"))]="DC"
  
  batch_lineage_celltype_NMF<-function(temp_celltype){
    Lineage_cells=rownames(Seuratobj@meta.data)[which(Seuratobj@meta.data$assign.curated  == temp_celltype)]
    Lineage_seurat=subset(Seuratobj,cells=Lineage_cells)
    NMF_list<<-list()
    batch_for_patient_celltype_NMF<-function(temp_sample){
      temp_cells=rownames(Lineage_seurat@meta.data)[which(Lineage_seurat@meta.data$Patient == temp_sample)]
      if(length(temp_cells) >100){
        temp_seurat=subset(Lineage_seurat,cells=temp_cells)
        temp_seurat=CreateSeuratObject(temp_seurat@assays$RNA@counts,meta.data=temp_seurat@meta.data)
        temp_seurat=NormalizeData(temp_seurat)
        temp_seurat=FindVariableFeatures(temp_seurat)
        temp_seurat=ScaleData(temp_seurat,do.center = FALSE)
        expmat=temp_seurat@assays$RNA@scale.data
        res <- nmf(expmat, 10,  seed = 10,method="snmf/r")
        saveRDS(res,file=paste0("/fs/home/hanya/Project/TME_Immune_difference/NatCancer_revision/ITH_NMF/",temp_dataset,"_",temp_celltype,"_",temp_sample,"_Prog10_NMF.rds"))
      }
    }
    temp_result=apply(as.matrix(unique(Lineage_seurat@meta.data$Patient)),1,batch_for_patient_celltype_NMF)
  }
  result=apply(as.matrix(unique(Seuratobj@meta.data$assign.curated)),1,batch_lineage_celltype_NMF)
}
result=apply(as.matrix(dataset),1,batch_process_nmf_program)


# for each cell lineage, summarized each NMF program using the top 50 genes based on NMF coefficients

NMF_file=list.files("./ITH_NMF",recursive = TRUE,full.names = TRUE)

Top_Feature<<-NULL
batch_process_NMF_file<-function(temp_file){
  temp_result=readRDS(temp_file)
  f <- extractFeatures(temp_result, 50L)
  f <- lapply(f, function(x) rownames(temp_result)[x])
  f <- do.call("rbind", f)
  f<-t(f)
  sample_name=gsub("_Prog10_NMF.rds","",unlist(strsplit(temp_file,"\\/"))[9])
  colnames(f)<-paste(sample_name,1:10)
  Top_Feature<<-cbind(Top_Feature,f)
}
result=apply(matrix(NMF_file),1,batch_process_NMF_file) 


#######Program similarity######
jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  #union = length(A) + length(B) - intersection
  return (intersection)
}
jaccrad_sim<<-NULL
batch_calculate_Jaccard_index<-function(temp_prog_index){
  reult=as.numeric(apply(matrix(1:ncol(Top_Feature)),1,function(y){ return(jaccard_similarity(Top_Feature[,temp_prog_index],Top_Feature[,y]))}))
  jaccrad_sim<<-rbind(jaccrad_sim,reult)
  #return(reult)
}
result<-apply(matrix(1:ncol(Top_Feature)),1,batch_calculate_Jaccard_index)

#######Program quality control::Process each program;#####
#non-redundant within the tumour
Sample_name=unique(unlist(lapply(strsplit(colnames(jaccrad_sim)," "),function(x) x[1])))
Valid_redundant_NMF_module<-function(temp_sam){
  index=grep(temp_sam,colnames(jaccrad_sim))
  temp_jaccard_sim=jaccrad_sim[index,index]
  valid_within_similarity<-function(temp_index){
    if(length(which(temp_jaccard_sim[temp_index,(temp_index+1):10] >10)) >0){
      return(colnames(temp_jaccard_sim)[which(temp_jaccard_sim[temp_index,(temp_index+1):10] >10)])
    }
  }
  result=apply(matrix(1:9),1,valid_within_similarity)
}
result=apply(matrix(Sample_name),1,Valid_redundant_NMF_module)
redundant_program=as.vector(unlist(result))
jaccrad_sim_subset=jaccrad_sim[match(setdiff(rownames(jaccrad_sim),redundant_program),rownames(jaccrad_sim)),match(setdiff(rownames(jaccrad_sim),redundant_program),colnames(jaccrad_sim))]
# Validation the robust of each program:

temp_result=apply(jaccrad_sim,1,function(x){length(which(x >= 10))})
temp_result_more10=names(temp_result)[which(temp_result > 1)]

jaccrad_sim_subset= jaccrad_sim_subset[match(temp_result_more10,rownames(jaccrad_sim_subset)),match(temp_result_more10,colnames(jaccrad_sim_subset))]
