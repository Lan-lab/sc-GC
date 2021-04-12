import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
sc.settings.set_figure_params( frameon=True,color_map='Spectral_r') 
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False

def data_process(data, gene_filter=True, min_genes=200, min_cells=10, min_counts=400, 
					max_genes=None,	max_genes_pct=99.5, mt_pct=0.3, mt_pct2=0.2,ncount2=1000,max_fraction=0.1,plot=True):
	#
	data.var_names_make_unique()
	#filter cells and genes 
	sc.pp.filter_cells(data,  min_counts=min_counts)
	sc.pp.filter_cells(data,  min_genes=min_genes)
	if gene_filter:
		sc.pp.filter_genes(data, min_cells=min_cells) 
	mito_genes = data.var_names.str.startswith('MT-')
	data.obs['percent_mito'] = np.sum(data[:, mito_genes].X, axis=1).A1 / np.sum(data.X, axis=1).A1
	data.obs['n_counts'] = data.X.sum(axis=1).A1
	if plot:
		sc.pl.violin(data, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
		sc.pl.scatter(data, x='n_counts', y='percent_mito')
		sc.pl.scatter(data, x='n_counts', y='n_genes')
	data = data[data.obs['percent_mito'] < mt_pct, :]
	data = data[  ~((data.obs['percent_mito'] > mt_pct2) & (data.obs['n_counts']<ncount2))  , :]
	if max_genes is None:
		max_genes=np.percentile(data.obs['n_genes'],max_genes_pct)
		print('n_genes at',max_genes_pct,'is',max_genes)
	data = data[data.obs['n_genes'] < max_genes, :]
	sc.pp.filter_cells(data,  min_counts=min_counts)
	sc.pp.filter_cells(data,  min_genes=min_genes)
	sc.pp.normalize_total(data,target_sum =None, exclude_highly_expressed =True, max_fraction = max_fraction)
	sc.pp.log1p(data, base=2)  #e
	data.raw=data 
	return data

# Load the raw data
exp_matrix = pd.read_csv('raw_exp_matrix.csv', index_col=0)
meta_matrix = pd.read_csv('exp_meta.csv', index_col=0)
sc_data = sc.Adata(exp_matrix)
sc_data.obs = meta_matrix

# Process the data
sc_data = data_process(sc_data, min_genes=200,min_counts=400, min_cells=10,mt_pct=0.3, mt_pct2=0.3,ncount2=1000, max_fraction =0.1)
sc.pp.highly_variable_genes(sc_data, n_top_genes = 6000, flavor ='seurat', batch_key='batch')
sc_scaled = sc.pp.scale(sc_data, max_value=10,copy=True)

# Data visualization
sc.tl.pca(sc_scaled,svd_solver='arpack', use_highly_variable = True)
sc.pp.neighbors(sc_scaled, n_pcs=10, n_neighbors=500)
sc.tl.umap(sc_scaled, init_pos=sc_scaled.obsm['X_pca'][:,:2] *0.0001)
sc.pl.umap(sc_scaled, color=['tissue'],size=1)
sc.pl.umap(sc_scaled, color=['patient'],size=1)
sc.pl.umap(sc_scaled, color=['cluster'],size=1)

# Differentially expressed genes
sc.tl.rank_genes_groups(sc_scaled, 'cluster', method='wilcoxon' , n_genes=500)
sc.tl.filter_rank_genes_groups(sc_scaled,max_out_group_fraction=0.6,min_in_group_fraction=0.1, min_fold_change = 2)
sc.pl.rank_genes_groups_heatmap(sc_scaled ,groupby='cluster',n_genes=5, swap_axes=True, figsize=(6,7),use_raw=False,cmap='RdBu_r')
sc.pl.rank_genes_groups_dotplot(sc_scaled ,groupby='cluster',n_genes=5, use_raw=False)














