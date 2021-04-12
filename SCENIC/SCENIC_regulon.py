from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
import os
sc.settings.set_figure_params( frameon=True,color_map='Spectral_r' )
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False

# Run time for demo data: several minutes.

DATA_FOLDER='./Data/'
exp_matrix = pd.read_csv(os.path.join(DATA_FOLDER,'DC_exp.csv'),index_col=0)
exp_meta = pd.read_csv(os.path.join(DATA_FOLDER,'DC_meta.txt'),sep='\t', index_col=0)

## Generate the regulons
Motifs_NAME = 'exp_matrix'
motifs = load_motifs(os.path.join(DATA_FOLDER, '{}.motifs.csv'.format(Motifs_NAME)))
regulons = df2regulons(motifs)

## Make a meta matrix for the regulon
reg_num = []
reg_target=[]
reg_tf=[]
for i in regulons:
	reg_tf.append(i.transcription_factor)
	reg_target.append(list(i.gene2weight.keys()))
	reg_num.append(len(list(i.gene2weight.keys())))

reg_meta = pd.DataFrame([reg_num,reg_target]  ).T
reg_meta.index=reg_tf
reg_meta.columns=['n_targets','targets']


## Calculate the AUCell scores
auc_mtx = aucell(exp_matrix, regulons, num_workers=8)
auc_mtx.columns = auc_mtx.columns.str[:-3]

## Analysis the result with scanpy
sc_auc_mtx = sc.AnnData(X=auc_mtx)
sc_auc_mtx.var_names =  [str(i) +'('+ str(j) + ')' for i, j in zip(auc_mtx.columns, reg_meta.n_targets)]
sc_auc_mtx.obs = exp_meta

## Differential analysis for AUCell scores
sc.tl.rank_genes_groups(sc_auc_mtx, 'leiden_cluster', method='wilcoxon',use_raw=True)
sc.tl.dendrogram(sc_auc_mtx, groupby='leiden_cluster')

sc.pl.rank_genes_groups_heatmap(sc_auc_mtx ,groupby='leiden_cluster',n_genes=5,swap_axes=True, figsize=(6,7),cmap='RdBu_r')
sc.pl.rank_genes_groups_dotplot(sc_auc_mtx ,groupby='leiden_cluster',n_genes=5)














