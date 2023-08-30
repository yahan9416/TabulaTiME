# @Author: Ya Han
# @Dat: 2020-08-06
# @Last Modified by: Ya
# @Last Modified time: 2020-12_10

# Note:
# The analysis pipeline starts from a count matrix, which can be MTX, HDF5 or plain text formatted, and then perform QC, clustering and cell-type annotation.
# The analysis pipeline generates a R script which can be modified, clustering and annotated UMAP plots and a Seurat object which can be readable and writable.
# Meta data is optional, not required.

# Acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/MAESTRO


# run the following command for detailed usage of tools.
MAESTRO mtx-to-h5 --help
MAESTRO merge-h5 --help
MAESTRO scrna-analysis --help

# the following commands can be used for format conversion
MAESTRO mtx-to-h5           # Convert 10X mtx format matrix to HDF5 format.
MAESTRO merge-h5            # Merge 10X HDF5 files
MAESTRO scrna-analysis      # scRNAseq data analysis pipeline


# Examples:
cd /home/hanya/Carcinogenesis/scRNAseq/P01_OSF/filtered_feature_bc_matrix
MAESTRO mtx-to-h5 --type Gene --matrix matrix.mtx --feature features.tsv --gene-column 2 --barcode barcodes.tsv --species GRCh38 --directory . --outprefix "P01_OSF"

MAESTRO merge-h5 --type Gene --h5 P01_OSF.h5 P01_Tumor.h5

cd  /home/hanya/Oral_cancer
MAESTRO scrna-analysis --format h5 --matrix ./Carcinogenesis_merge_delP12OLP_addP13_gene_count.h5 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix Carcinogenesis_HNSCC

