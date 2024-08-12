# @Author: Ya Han
# @Dat: 2024-08-12
# @Last Modified by: Ya



####Acivate conda environment
source activate ./miniconda3/envs/MAESTRO


####Quality control
# Convert 10X mtx format matrix to HDF5 format.
MAESTRO mtx-to-h5 --type Gene --matrix matrix.mtx --feature features.tsv --gene-column 2 --barcode barcodes.tsv --species GRCh38 --directory . --outprefix "Sample1"
# Merge 10X HDF5 files
MAESTRO merge-h5 --type Gene --h5 sample1.h5 sample2.h5
####scRNAseq data analysis pipeline
#The output is Seurat object
MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_GSE1.h5 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix TabulaTIME_GSE1


##Doublets 
#Doublets prediction for each samples
Python Doublet_scrublet.py Need_process_file_path
#Doublets removal
R Doublet_remove_from_seuratobj.R Scrublet_result_path Seurat_object_path

### Data Pre-processing
####Malignant cell identification
#Mannually selected one Copy number variation prediction methods
#CopyKat
R Copykat.R RowCount_matrix
# inferCNV 
R inferCNV.R RowCount_matrix inferCNV_gene_order_file.txt


####Batch effect
#Evaluating the batch effect
R Entropy_evaluated_batcheffect.R Seurat_object_path
#Corrected batch effect
R Correct_batch_effect_by_CCA.R Seurat_object_path
