# @Author: Ya Han
# @Dat: 2024-08-12
# @Last Modified by: Ya

####Acivate conda environment
source activate ./miniconda3/envs/MAESTRO

# ST Data collection and pre-processing
#QC
R Process_ST_data.R ST_Feature_bc_matrix_path

#Malignant cell identification
R ST_copykat.R ST_seurat

#Spot deconvolutaion 
STRIDE deconvolve --sc-count scRNA_count_matrix.txt \
                  --sc-celltype scRNA_meta_info.txt \
                  --st-count  filtered_feature_bc_matrix.h5 \
                  --outdir Result/STRIDE --outprefix Sample --normalize

