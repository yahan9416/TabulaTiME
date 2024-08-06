import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import re

pattern = "count_matrix"
os.chdir('./Doublets_Finder/Scrublet')
for file in os.listdir('./Doublets_Finder'):
	if(re.search(pattern,file)):
		counts_matrix = scipy.io.mmread('./Doublets_Finder/'+file).T.tocsc()
		scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
		doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
		np.savetxt("dbl_" + file + ".txt", scrub.predicted_doublets_, fmt='%d')