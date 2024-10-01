#' Deconvoluting the bulk RNA-seq expression matrix using EPIC.
#' @param exp_path Path of expression matrix
#' @author Ya Han


args <- commandArgs(trailingOnly = TRUE)
exp_path <- args[1]


library(EPIC)
microarray=read.table(exp_path,header = TRUE,sep="\t",row.names = 1)
res3 <- EPIC(bulk=microarray, reference=TRef)
names(res3)
head(res3$cellFractions)
head(res3$cellFractions$CD8_Tcells)
