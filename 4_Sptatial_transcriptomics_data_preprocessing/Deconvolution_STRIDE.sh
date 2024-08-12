STRIDE deconvolve --sc-count scRNA_count_matrix.txt \
                  --sc-celltype scRNA_meta_info.txt \
                  --st-count  filtered_feature_bc_matrix.h5 \
                  --outdir Result/STRIDE --outprefix Sample --normalize

