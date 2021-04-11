import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import sparse

scv.set_figure_params(color_map='RdBu_r', frameon=True)
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False


#load the data
DC_spliced = sparse.load_npz('./Data/DC_spliced.npz') 
DC_unspliced = sparse.load_npz('./Data/DC_unspliced.npz') 
DC_normed = pd.read_csv('./Data/DC_normed.txt',sep='\t', index_col=0)
DC_meta = pd.read_csv('./Data/DC_meta.txt',sep='\t', index_col=0)

sc_DC = sc.AnnData(DC_normed)
sc_DC.obs = DC_meta
sc_DC.layers['spliced'] = DC_spliced
sc_DC.layers['unspliced'] = DC_unspliced
scv.pp.filter_and_normalize(sc_DC, min_cells =3, min_cells_u=3,min_shared_cells=3, n_top_genes=None, flavor='seurat', log=True, copy=False)


sc.tl.pca(sc_DC,svd_solver='arpack')
sc.pp.neighbors(sc_DC,n_pcs=20, n_neighbors=30)
sc.tl.umap(sc_DC,init_pos=sc_DC.obsm['X_pca'][:,:2] *0.0001)
sc.pl.umap(sc_DC, color=['leiden_cluster'],size=10)



# Calculate the RNA velocity
scv.pp.moments(sc_DC)
scv.tl.recover_dynamics(sc_DC,n_jobs=8)
scv.tl.velocity(sc_DC, mode='dynamical')
scv.tl.velocity_graph(sc_DC)

scv.pl.velocity_embedding(sc_DC, basis='umap', color = 'tissue', arrow_size = 3, arrow_length=3)
scv.pl.velocity_embedding_stream(sc_DC, basis='umap', color = 'leiden_cluster',legend_fontsize=5)


# Calculate the pseudotime (called latent time here)
scv.pl.velocity_embedding(sc_DC, basis='umap', color = 'leiden_cluster', arrow_size = 3, arrow_length=1)

scv.pl.velocity_embedding_stream(sc_DC, basis='umap', color = 'leiden_cluster',legend_fontsize=12)

scv.tl.recover_latent_time(sc_DC)

scv.pl.scatter(sc_DC, color='latent_time',basis='umap',cmap='YlGnBu')




