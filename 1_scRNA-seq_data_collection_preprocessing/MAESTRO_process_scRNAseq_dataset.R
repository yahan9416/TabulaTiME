#' Porcess the cell expression by MAESTRO
#' @param Expmat_path Path of expression matrix
#' @param ncores Number of cores used 
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Expmat_path<-args[1]
ncores<-args[2]

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 10*1024^3)

# read data
expr = Read10X_h5(Expmat_path)

# choose optimal pc and resolution based on cell number
cells <- ncol(expr)
if (cells <= 5000) {
  dims.use <- 1:15; cluster.res <- 0.6
} else if(cells <= 10000) {
  dims.use <- 1:20; cluster.res <- 1
} else if(cells <= 40000) {
  dims.use <- 1:30; cluster.res <- 1
} else if(cells <= 80000) {
  dims.use <- 1:40; cluster.res <- 1
} else if(cells <= 150000) {
  dims.use <- 1:50; cluster.res <- 1
} else {
  dims.use <- 1:75; cluster.res <- 1
}

# choose npc
npc <- ifelse(max(dims.use) < 50, 50, 
              ifelse(max(dims.use) < 75, 75, 100))

# clustering
RNA.res = RNARunSeurat(inputMat = expr, 
                       project = "TabulaTIME_GSE1", 
                       min.c = 10,
                       min.g = 500,
                       runpca.agrs = list(npcs = npc),
                       dims.use = dims.use,
                       variable.genes = 2000, 
                       organism = "GRCh38",
                       cluster.res = 0.6,
                       genes.test.use = "presto",
                       only.pos = TRUE,
                       genes.cutoff = 1e-05)

# cell-type annotation
RNA.res$RNA = RNAAnnotateCelltype(RNA = RNA.res$RNA, 
                                  genes = RNA.res$genes,
                                  signatures = "human.immune.CIBERSORT",
                                  min.score = 0.1)

# save object
saveRDS(RNA.res, "TabulaTIME_GSE1_seurat.rds")
