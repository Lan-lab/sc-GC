import numpy as np
import pandas as pd
import seaborn as sns

## Excepted results are TCR-clonotype sharing matrices for CD4 T cells and CD8 T cells, respectively.

## Functions
def TCR_process(T_meta, TCR_col, cluster_col, cluster_names):
	clone_type_table = T_meta[TCR_col].value_counts()
	clone_type_names = clone_type_table.index
	for i in cluster_names:
		obs_one_cluster = T_meta[T_meta[cluster_col]==i]
		one_cluster_table = obs_one_cluster[TCR_col].value_counts()
		clone_type_table = pd.concat([clone_type_table,one_cluster_table],axis=1)
	clone_type_table = clone_type_table.fillna(0)
	clone_type_table.columns = ['total']+cluster_names 
	clone_type_table_part = clone_type_table.reindex(columns=cluster_names)
	df_filtered = clone_type_table_part[clone_type_table_part.max(axis=1)>=3].copy(deep=True)
	# get peak position
	peak_position=[]
	for i in range(len(df_filtered)):
		peak_position.append( df_filtered.iloc[i].idxmax() )
	df_filtered['peak_position'] = list(pd.Series(peak_position))
	return df_filtered
	
def TCR_sharing(clone_type_table_peak_position, sharing_threshold=3):
	clone_type_distribution = clone_type_table_peak_position.peak_position.value_counts()
	clusters = clone_type_table_peak_position.columns[:-1]
	TCR_sharing_matrix=[]
	for i in clusters:
		if i in clone_type_distribution.index:
			TCR_sharing_vector=[]
			clone_type_table_sub = clone_type_table_peak_position[clone_type_table_peak_position.peak_position==i]
			for j in clusters:
				shared_clone = np.count_nonzero(clone_type_table_sub[j]>=sharing_threshold)
				TCR_sharing_vector.append(  shared_clone/clone_type_distribution[i])
			TCR_sharing_matrix.append(TCR_sharing_vector)
		else:
			TCR_sharing_matrix.append([0]*len(clusters))
	TCR_sharing_df = pd.DataFrame(TCR_sharing_matrix,index=clusters,columns=clusters)
	for i in range(len(clusters)):
		TCR_sharing_df.iloc[i,i] = np.nan
	return TCR_sharing_df


### Two subsets of T cells.
CD8_cluster_names = ['Proliferative_T_cells',
 'CD8_C1_LEF1',
 'CD8_C2_CX3CR1',
 'CD8_C3_GZMH',
 'CD8_C4_GZMK',
 'CD8_C5_TOB1',
 'CD8_C6_GNLY',
 'CD8_C7_CD160',
 'CD8_C8_IL17A',
 'CD8_C9_HAVCR2',
 'CD8_C10_SLC4A10']
CD4_cluster_names = ['Proliferative_T_cells',
 'CD4_C1_CCR7',
 'CD4_C2_LTB',
 'CD4_C3_SLC2A3',
 'CD4_C4_CD69',
 'CD4_C5_CXCL13',
 'CD4_C6_IL17A',
 'Treg_C1_SELL',
 'Treg_C2_LAG3',
 'Treg_C3_CTLA4']

T_meta = pd.read_csv('./Data/TCR_meta.txt',sep='\t')

CD8_clone_type_table = TCR_process(T_meta, TCR_col='TCR_merge', cluster_col = 'leiden_cluster', cluster_names = CD8_cluster_names )
CD4_clone_type_table = TCR_process(T_meta, TCR_col='TCR_merge', cluster_col = 'leiden_cluster', cluster_names = CD4_cluster_names )

#TCR-clonotype sharing matrix of CD8 T cells
CD8_TCR_sharing_df = TCR_sharing(CD8_clone_type_table, sharing_threshold=3)
#TCR-clonotype sharing matrix of CD4 T cells
CD4_TCR_sharing_df = TCR_sharing(CD4_clone_type_table, sharing_threshold=3)

#Plot
ax = sns.clustermap(CD8_TCR_sharing_df, annot=True, fmt=".2f",cmap='RdBu_r',figsize=(5,5),row_cluster=False , col_cluster=False)
ax.fig.autofmt_xdate(rotation=45, ha='right')

ax = sns.clustermap(CD4_TCR_sharing_df, annot=True, fmt=".2f",cmap='RdBu_r',figsize=(5,5),row_cluster=False , col_cluster=False)
ax.fig.autofmt_xdate(rotation=45, ha='right')


