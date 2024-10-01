#' Deconvolution the ST data by RCTD
#' @param Sc_Seurat_path Path of single cell seurat object
#' @param ST_Seurat_path Path of ST seurat object
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Sc_Seurat_path<-args[1]
ST_Seurat_path<-args[2]

library(Seurat)
library(MAESTRO)
library(spacexr)

#Reference building
##############

P13_Seurat=readRDS(Sc_Seurat_path)
table(P13_Seurat@meta.data$Source)
DefaultAssay(P13_Seurat)<-"RNA"
Nor_cells_index=rownames(P13_Seurat@meta.data)[which(P13_Seurat@meta.data$Source == "Normal")]
P13_Seurat_nor=subset(P13_Seurat,cells=Nor_cells_index)


less25_cells=which(table(P13_Seurat@meta.data$curated_anno_gene) <25)
doulet=rownames(P13_Seurat@meta.data)[which(P13_Seurat@meta.data$curated_anno_gene %in% names(less25_cells))]
P13_Seurat=subset(P13_Seurat,cells=doulet,invert=TRUE)
#SetAssayData(P13_Seurat)<-"counts"
count.data <- GetAssayData(object = P13_Seurat[["RNA"]], slot = "counts")
Meta_Info=P13_Seurat@meta.data
cell_types=Meta_Info$curated_anno_gene
names(cell_types) <- rownames(Meta_Info)
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- Meta_Info$nCount_RNA
names(nUMI) <- rownames(Meta_Info)
reference <- Reference(count.data, cell_types, nUMI)
#saveRDS(reference,'ref.rds')

#######ST process#########

P13_ST_seurat=readRDS(ST_Seurat_path)
DefaultAssay(P13_ST_seurat)<-"Spatial"
coords <- GetTissueCoordinates(P13_ST_seurat, scale = NULL)
P13_ST_counts <- as.matrix(GetAssayData(P13_ST_seurat, assay = "Spatial", slot = "counts"))
nUMI <- P13_ST_seurat@meta.data$nCount_Spatial
names(nUMI)<-rownames(P13_ST_seurat@meta.data)
### Create SpatialRNA object
puck <- SpatialRNA(coords, P13_ST_counts, nUMI)



library(spacexr)
myRCTD <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE) # here puck is the SpatialRNA object, and reference is the Reference object.
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
saveRDS(myRCTD,file="RCTD_deconvolution_result.rds")


#####Process RCTD result##########
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
#resultsdir <- 'RCTD_RefAdd_Plots'
resultsdir <- 'RCTD_Ref_CorrEpi_Plots'
#resultsdir <- 'RCTD_RefDel_Plots' ## you may change this to a more accessible 
#directory on your computer.
dir.create(resultsdir)
# make the plots 
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 
