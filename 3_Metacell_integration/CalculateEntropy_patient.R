#' Calculating the Entropy score to evaluate the effect of batch effect
#' @param data_dirs The path of process scRNA-seq Seurat object
#' @author Ya Han

library(Seurat)
library(MAESTRO)
library(rjson)
library(dbscan)


args <- commandArgs(trailingOnly = TRUE)
data_dirs<-args[1]

CalEntropy30 <- function(rds,Vbatch){
  # k <- ceiling(dim(rds)[2]/100)
  k = 30
  knn <- kNN(rds@reductions[["umap"]]@cell.embeddings,k = k)
  neighbor_matrix <- knn$id
  #fetch a vector of batch conditions: ["patient1","patient1","patient2"...]
  #Vbatch <- rds@meta.data[,BatchName]
  #batch compositions: ["patient1","patient2","patient3"...]
  batch <- unique(Vbatch)
  Entropies = sapply(1:ncol(rds), function(i){
    sum=0
    cell_neighbors <- Vbatch[neighbor_matrix[i,]]
    batch_pct = sapply(1:length(batch), function(j){
      Nbatch <- length(which(cell_neighbors==batch[j]))
      if (Nbatch == 0){
        return(0)
      } else {
        return((Nbatch/k)*log2(Nbatch/k))
      }
    })
    return(-sum(batch_pct))
  })
  # return entropies(a vector of number): [1,0.5,0.3 ...]
  return(Entropies)
}

CalEntropy <- function(rds,Vbatch){
  k <- ceiling(dim(rds)[2]/100)
  # k = 30
  knn <- kNN(rds@reductions[["umap"]]@cell.embeddings,k = k)
  neighbor_matrix <- knn$id
  #fetch a vector of batch conditions: ["patient1","patient1","patient2"...]
  #Vbatch <- rds@meta.data[,BatchName]
  #batch compositions: ["patient1","patient2","patient3"...]
  batch <- unique(Vbatch)
  Entropies = sapply(1:ncol(rds), function(i){
    sum=0
    cell_neighbors <- Vbatch[neighbor_matrix[i,]]
    batch_pct = sapply(1:length(batch), function(j){
      Nbatch <- length(which(cell_neighbors==batch[j]))
      if (Nbatch == 0){
        return(0)
      } else {
        return((Nbatch/k)*log2(Nbatch/k))
      }
    })
    return(-sum(batch_pct))
  })
  # return entropies(a vector of number): [1,0.5,0.3 ...]
  return(Entropies)
}

entropy_dir = "/home1/wangchenfei/Project/TIRA/Data_Entropy/Patient_CCA"

patient_file = "Batch_patient.txt"
CCA_patient_project = readLines(patient_file)
for (data_dir in data_dirs) {
  setwd(data_dir)
  dirs = list.files()
  for (dir in dirs) {
    objects = list.files(dir, "_res.rds")
    if (length(which(grepl("CCA", objects))) > 0){
      objects = objects[grepl("CCA", objects)]
    }
    for (i in objects){
      project = gsub("_CCA", "", i)
      project.name = gsub("_res.rds", "", project)
      if (project.name %in% CCA_patient_project) {
        SeuratObj = readRDS(file.path(dir,i))
        project.name = SeuratObj$RNA@project.name
        
        project.name = gsub("_CCA", "", project.name)
        metafile = file.path(dir, paste0(project.name, "_umap.json"))
        meta_list = fromJSON(file = metafile)
        meta_df = as.data.frame(meta_list)
        rownames(meta_df) = meta_df$Cell
        batch = "Patient"
        Vbatch = meta_df[colnames(SeuratObj$RNA), batch]
        entropies30 = CalEntropy30(SeuratObj$RNA, Vbatch)
        entropies = CalEntropy(SeuratObj$RNA, Vbatch)
        res_df = data.frame(Cell = colnames(SeuratObj$RNA), entropy30 = entropies30, entropy = entropies, stringsAsFactors = FALSE)
        write.table(res_df, file = file.path(entropy_dir, paste0(project.name, "_", batch, "_", length(unique(Vbatch)), "_CCA_cell_entropy.txt")),
                    sep = "\t", quote = FALSE)
        message(i)
      }
    }
  }
}
