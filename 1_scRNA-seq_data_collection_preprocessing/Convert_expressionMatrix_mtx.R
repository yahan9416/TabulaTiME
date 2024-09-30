#' Writing the expression matrix into MTX format.
#' @param Minicluster_TPM Expression matrix
#' @param dir Path of file store
#' @param file_name names of output MTX file
#' @author Ya Han
#' 
#The order of args: Expression matrix, Direction of matrix writein, file_name
args <- commandArgs(trailingOnly = TRUE)
Minicluster_TPM <- args[1]
dir <- args[2]
file_name <- args[3]

Minicluster_TPM <- as(Minicluster_TPM, "sparseMatrix")
mtx_file = file.path(dir, paste0(file_name, "_count_matrix.mtx"))
gene_file = file.path(dir, paste0(file_name, "_count_genes.tsv"))
barcode_file = file.path(dir, paste0(file_name, "_count_barcodes.tsv"))
writeMM(Minicluster_TPM, mtx_file)
write.table(rownames(Minicluster_TPM), gene_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(colnames(Minicluster_TPM), barcode_file, col.names = FALSE, row.names = F, quote = FALSE)
