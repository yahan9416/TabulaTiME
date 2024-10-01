## Run_Metacell_integration.sh

### Step 1. MetaCell integration
    #Merge MetaCell expression matrix
    MAESTRO merge-h5 --type Gene --h5 Datase1_Metacell_exp.h5 Datase2_Metacell_exp.h5

    #scRNA-seq data pre-processing 
    MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_megred_metacell_exp.h5 count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 directory . --outprefix TabulaTIME_megred_metacell_exp

### Step 2. Batch effect correction   
    #Correct batch effect by CCA
    R CCA_corrected_batch_effect_for_all_Metacell.R Seurat_obj_path Batch_info Nfeature N_dimused ncores

### Step 3. Evaluating the performance of batch effect correction
    #Local Inverse Simpson's Index (LISI)
    #Calculating the LISI score to evaluate the performance of batch effect correction
    R LISI_quantitatively_evaluate_integration.R CCA_seurat_obj_path
    #Entropy
    #Calculating the Entropy score to evaluate the performance of batch effect correction.
    R Entropy_evaluated_batcheffect.R CCA_seurat_obj_path
    #Adjusted Rand Index (ARI) was used to evaluate the performance of batch effect correction.
    R ARI_adjusted_rank_index.R CCA_seurat_obj_path
    #Calculating the average silhouette score (ASW) for each cell to evaluate the optimal cluster resolution.
    R Silhouette_ASW.R CCA_seurat_obj_path

### Step 4. Identifying the optimal clustering resolution
    #Calculating the average silhouette score for each cell to evaluate the optimal cluster resolution.
    R Silhouette_ASW.R CCA_seurat_obj_path
    #Clustree were used to select the optimal cluster resolution 
    R Clustree.R CCA_seurat_obj_path Marker_genes

### Step 5. Metacell annotation
    #Calculate the Rogue score to estimate the purity of each cluster or cell type.
    R Rogue.R CCA_seurat_obj
