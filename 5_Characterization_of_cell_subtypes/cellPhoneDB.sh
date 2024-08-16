# @Author: Ya Han
# @Dat: 2021-07-01
# @Last Modified by: Ya


# Note:
# The analysis pipeline starts from the expression matrix of cells and the cell type information of each cell type


# Acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/cpdb


# run the following command for detailed usage of tools.
cellphonedb --help
cellphonedb method --help
cellphonedb method statistical_analysis --help



# Examples:
cd /fs/home/hanya/Project/Carcinogenesis/scRNAseq/Cell_Cell_interacton/CD8_Exhausted_interaction/After_RNAvelocity_reHeal
nohup cellphonedb method statistical_analysis ./OX40L_CD8Havcr2_minicluster_meta.txt ./OX40L_CD8Havcr2_celltype_minicluster_exp_count.txt --threads=40 --counts-data gene_name --project-name=OX40L_CD8Havcr2_celltype_inter >log.txt & 
      
