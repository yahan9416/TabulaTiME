#Deconvolution by RCTD
library(Seurat)
library(MAESTRO)
library(spacexr)

#Reference building
##############
setwd("/Storage/hanya/Spatial_process/")
P13_Seurat=readRDS("P13_Seurat_subcell_add_less25cells.rds")
table(P13_Seurat@meta.data$Source)
meta_data=readRDS("P13_Seurat_subcell_add_less25cells_meta.rds")
P13_Seurat@meta.data=meta_data
DefaultAssay(P13_Seurat)<-"RNA"
Nor_cells_index=rownames(P13_Seurat@meta.data)[which(P13_Seurat@meta.data$Source == "Normal")]
P13_Seurat_nor=subset(P13_Seurat,cells=Nor_cells_index)
predict=read.table("/Storage/hanya/Spatial_process/test_copykat_prediction.txt",header = TRUE,sep="\t")
epi_p13=intersect(predict[which(predict[,2] == "diploid"),1],Nor_cells_index)
P13_Seurat_nor@meta.data$curated_anno_gene[which(rownames(P13_Seurat_nor@meta.data) %in% epi_p13)]="Epithelial"
sort(table(P13_Seurat_nor@meta.data$curated_anno_gene))


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
setwd("/fs/home/hanya/Project/Carcinogenesis/Spatial_processed/P13")
setwd("/Storage/hanya/Spatial_process")
P13_ST_seurat=readRDS("Oral_ST_P13Normal_seurat.rds")
DefaultAssay(P13_ST_seurat)<-"Spatial"
coords <- GetTissueCoordinates(P13_ST_seurat, scale = NULL)
P13_ST_counts <- as.matrix(GetAssayData(P13_ST_seurat, assay = "Spatial", slot = "counts"))
nUMI <- P13_ST_seurat@meta.data$nCount_Spatial
names(nUMI)<-rownames(P13_ST_seurat@meta.data)
### Create SpatialRNA object
puck <- SpatialRNA(coords, P13_ST_counts, nUMI)


setwd("/Storage/hanya/Spatial_process/RCTD")
library(spacexr)
myRCTD <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE) # here puck is the SpatialRNA object, and reference is the Reference object.
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
saveRDS(myRCTD,file="Oral_ST_P13Normal_Ref_CorrectMali_result_full.rds")
saveRDS(myRCTD,file="Oral_ST_P13Normal_RefAdd_RCTD_result_full.rds")
saveRDS(myRCTD,file="Oral_ST_P13Normal_RefDel_RCTD_result_full.rds")

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
