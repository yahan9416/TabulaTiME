# @Author: Ya Han
# @Dat: 2024-08-12
# @Last Modified by: Ya

####Acivate conda environment
source activate ./miniconda3/envs/MAESTRO

####Metacell Integration
MAESTRO merge-h5 --type Gene --h5 Datase1_Metacell_exp.h5 Datase2_Metacell_exp.h5
MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_megred_metacell_exp.h5 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix TabulaTIME_megred_metacell_exp

####Batch effect correction
R CCA_corrected_batch_effect_for_all_Metacell.R Integrated_seurat_obj

####Evalutation the optimal batch effect correction method
#LISI
R LISI_quantitatively_evaluate_integration.R CCA_seurat_obj
#Entropy
R Entropy_evaluated_batcheffect.R CCA_seurat_obj
#ASW
R Silhouette_ASW.R CCA_seurat_obj
#ARI 
R ARI_adjusted_rank_index.R CCA_seurat_obj