#' Running CopyKat to predict malignant cells.
#' @param rowcount_matrix Expression matrix
#' @param num_cores Number of core used
#' @author Ya Han

library(copykat)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
rowcount_matrix <- args[1]
num_cores <-args[2]

expmat=Read10X_h5(rowcount_matrix)#the path of raw count expression
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=num_cores)
predict=copykat.test$prediction[,2]