# @Author: Ya Han
# @Dat: 2024-07-01
# @Last Modified by: Ya


# Note:
# The analysis pipeline starts from the expression matrix of cells and the cell type information of each cell type


# Acivate conda environment

# run the following command for detailed usage of tools.
cellphonedb --help
cellphonedb method --help
cellphonedb method statistical_analysis --help



# Examples:
nohup cellphonedb method statistical_analysis ./meta_infp.txt ./exp_count.txt --threads=40 --counts-data gene_name --project-name=celltype_inter >log.txt & 
      
