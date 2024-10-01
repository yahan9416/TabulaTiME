#' Infer cell communication using CellChat
#' @param Seurat_obj_path Path of seurat object
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj_path<-args[1]


library(Seurat) # please update to Seurat V4
library(tidyverse)
library(MAESTRO)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


CD8Treg_seu=readRDS(Seurat_obj_path)

table(CD8Treg_seu@meta.data$curated_anno_gene)
normal_exp=GetAssayData(CD8Treg_seu,slot="data")
meta_infp=CD8Treg_seu@meta.data

cellchat <- createCellChat(object = normal_exp, meta = meta_infp, group.by = "curated_anno_gene")
rm(CD8Treg_seu)
gc()
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 5)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("Cellchat_Treg_CD8T.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("Cellchat_Treg_CD8T_seperate_each_celltype.pdf")
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show <- c("MIF")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,3)
#pathways.show.all=c("CXCL","CD45","CD80","CD86","OX40","CD40","NOTCH","WNT","MHC-I","IL4","IL16","CX3C")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
