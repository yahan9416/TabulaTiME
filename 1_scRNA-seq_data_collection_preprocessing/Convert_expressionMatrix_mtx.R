Minicluster_TPM <- as(Minicluster_TPM, "sparseMatrix")
dir="/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Integrated_analysis/Meger_minicluster_expression"
mtx_file = file.path(dir, paste0("Minicluster_Merge", "_count_matrix.mtx"))
gene_file = file.path(dir, paste0("Minicluster_Merge", "_count_genes.tsv"))
barcode_file = file.path(dir, paste0("Minicluster_Merge", "_count_barcodes.tsv"))
writeMM(Minicluster_TPM, mtx_file)
write.table(rownames(Minicluster_TPM), gene_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(colnames(Minicluster_TPM), barcode_file, col.names = FALSE, row.names = F, quote = FALSE)
