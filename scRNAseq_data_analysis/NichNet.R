#NichNET
library(nichenetr)
library(Seurat) 
library(tidyverse)
library(MAESTRO)

#Need prepare object
Fibro_Seurat
Fibro_DEG
intersted_celltype

#Need download from Nichnet
ligand_target_matrix=readRDS("ligand_target_matrix.rds")
lr_network=readRDS("lr_network.rds")
weighted_networks=readRDS("weighted_networks.rds")

DefaultAssay(Fibro_Seurat)<-"RNA"
expressed_genes_sender=rownames(Fibro_Seurat)
geneset_oi=rownames(Fibro_DEG)[which(Fibro_DEG$cluster == intersted_celltype)]
background_expressed_genes=expressed_genes_sender %>% .[. %in% rownames(ligand_target_matrix)]



ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_sender)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(100, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
write.table(ligand_activities,file="/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Lineage_seperate_analysis/Fibroblast/NichNet/NichNet_FibroCTHRC1_ligand_activities.txt",col.names = TRUE,sep="\t",quote = FALSE)
head(best_upstream_ligands)

# Infer target genes of top-ranked ligands and visualize in a heatmap

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
active_ligand_target_links_df=active_ligand_target_links_df[which(active_ligand_target_links_df$target %in% rownames(active_ligand_target_links)),]
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
pdf("Predicted_Upgene_ligand.pdf",width = 12,height = 9)
p_ligand_target_network
dev.off()

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

pdf("Celltype_ligand_receptor_network.pdf",width = 25,height = 9)
p_ligand_receptor_network
dev.off()


ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
pdf("Celltype_ligand_pearson.pdf",width = 3,height = 9)
p_ligand_pearson
dev.off()
