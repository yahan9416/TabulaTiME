# @Author: Ya Han
# @Dat: 2024-08-12
# @Last Modified by: Ya


####Acivate conda environment
source activate ./miniconda3/envs/MAESTRO


####Generating metacell within each datasets 
R Generate_metacell_for_each_dataset.R Seurat_object_path Cel_number_within_each_metacell


####Determining the optimal number of cells per metacell
#Gene coverage: Metacell_info_path is the output file of Generate_metacell_for_each_dataset.R
R Evaluated_gene_coverage.R Metacell_info_path
#Within metacell variation (Gini Index)
R GINI_index_within_metacell_variation.R Metacell_info_path Seurat_object_path


####Evaluating the performance of metacell
#Integrated metacell and generated Seurat object
MAESTRO merge-h5 --type Gene --h5 Datase1_Metacell_exp.h5 Datase2_Metacell_exp.h5
MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_megred_metacell_exp.h5 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix TabulaTIME_megred_metacell_exp
#LISI
R LISI_quantitatively_evaluate_integration.R Integrated_seurat_obj
#Biological program signature scores
R Biological_program_generated_NMF.R Seurat_object_path
