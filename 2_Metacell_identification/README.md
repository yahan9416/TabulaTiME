## Run_ MetaCell_identification.sh

### Step 1. Constructing MetaCells from scRNA-seq data
    #Generating metacell within each datasets 
    R Generate_metacell_for_each_dataset.R Seurat_object_path Cel_number_within_each_metacell
    
### Step 2. Determining the optimal number of cells per MetaCel
    #Gene coverage: Metacell_info_path is the output file of Generate_metacell_for_each_dataset.R
    R Evaluated_gene_coverage.R Metacell_info_path
    
    #Within metacell variation (Gini Index) 
    #Calculating the GINI index to evaluate the optimal number of cells within each metacell.
    R GINI_index_within_metacell_variation.R Metacell_info_path Seurat_object_path

### Step 3. Evaluating the performance of MetaCell
    #Local inverse Simpson’s Index (LISI)
    #The Local Inverse Simpson Index (LISI) score was calculated to evaluate the performance of batch effect correction.
    R LISI_quantitatively_evaluate_integration.R Integrated_seurat_obj
