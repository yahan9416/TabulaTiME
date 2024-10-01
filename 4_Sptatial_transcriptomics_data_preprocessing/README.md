## Run_ST_data_preprocesing.sh

### Step 1. Data collection

From the GEO database, we collect the ST data of 62 patients of 6 cancer types and download the gene-spot matrixes, cell coordinates, and corresponding tissue image

### Step 2. Quality control
    #UMI count per spot (>1000); Gene number per spot (>250); Precent MT gene per spot (<25%)
    R Process_ST_data.R ST_Feature_bc_matrix_path ncores
    
### Step 3. Data Pre-processing
    #Malignant cell identification: Copy number variation
    R ST_copykat.R ST_seurat_path ncores
    
    # Spot deconvolutaion
    STRIDE deconvolve --sc-count scRNA_count_matrix.txt --sc-celltype scRNA_meta_info.txt --st-count  filtered_feature_bc_matrix.h5  --outdir Result/STRIDE --outprefix Sample --normalize

### Step 4. Spatial localization analysis
    #Evaluating the cell type colocalozation in ST dataset 
    R Celltype_colocalization.R Celltype1_DEG_gene_list_path Celltype2_DEG_gene_list_path Deconvolution_result_path Seurat_obj_path
