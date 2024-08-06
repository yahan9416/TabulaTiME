cd /fs/home/hanya/Project/Carcinogenesis/Spatial_processed/P13/STRIDE/

STRIDE deconvolve --sc-count P13_Seurat_Nor_count_matrix.txt \
                  --sc-celltype P13_Seurat_Nor_meta_info.txt \
                  --st-count /fs/home/hanya/Project/Carcinogenesis/Spatial_processed/Data/P13NormalOLK_B1/filtered_feature_bc_matrix.h5 \
                  Deconvlution--outdir Result/STRIDE --outprefix Oral_P13Nor --normalize