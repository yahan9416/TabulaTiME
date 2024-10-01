## Run_Metacell_integration.sh

### Step 1. MetaCell integration
    #Merge MetaCell expression matrix
    MAESTRO merge-h5 --type Gene --h5 Datase1_Metacell_exp.h5 Datase2_Metacell_exp.h5

    #scRNA-seq data pre-processing 
    MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_megred_metacell_exp.h5 count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 directory . --outprefix TabulaTIME_megred_metacell_exp

### Step 2. Batch effect correction   
    #CCA
    R CCA_corrected_batch_effect_for_all_Metacell.R Integrated_seurat_obj_path

### Step 3. Evaluating the performance of batch effect correction
    #Local Inverse Simpson's Index (LISI)
    R LISI_quantitatively_evaluate_integration.R CCA_seurat_obj_path
    #Entropy
    R Entropy_evaluated_batcheffect.R CCA_seurat_obj_path
    #Adjusted rand index (ARI)
    R ARI_adjusted_rank_index.R CCA_seurat_obj_path
    #Average Silhouette Width (ASW)
    R Silhouette_ASW.R CCA_seurat_obj_path

### Step 4. Identifying the optimal clustering resolution
    #Silhouette score
    R Silhouette_ASW.R CCA_seurat_obj_path
    #Clustree
    R Clustree.R CCA_seurat_obj_path Marker_genes

### Step 5. Metacell annotation
    #ROGUE
    R Rogue.R CCA_seurat_obj
