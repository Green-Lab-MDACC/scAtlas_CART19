import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=300, facecolor='white')
# sc.settings.figdir = './figures/' # default


# reload
adata = sc.read('ana.h5ad')

# markers
markers = {'CD3': ['CD3D', 'CD3E', 'CD3G', 'CD247'],
           'CD8': ['CD8A', 'CD8B'],
           'CD4': ['CD4'],
           'IAC': ['IL1B', 'SIRPA', 'CXCL8', 'LILRB4']}

def knn_smooth_expr(adata, markers):
    """
    This function imputes expressions of markers using knn graph.

    Parameters:
    - adata: AnnData object containing log-normalized expressions stored in adata.raw and embeddings stored in adata.obsm['X_pca'].
    - markers: Dict with keys for marker names and values for gene lists.

    Returns:
    AnnData object containing knn smoothed expressions stored in adata.obs.

    Example:
    ```
    markers = {'CD8': ['CD8A', 'CD8B'], 'CD4': ['CD4']}
    adata = knn_smooth_expr(adata, markers)
    ```
    """

    # hack sc.pp.neighbors that results in
    # connectivities: [0, 1] values calculated from umap.fuzzy_simplicial_set using knn_indices and knn_dists;
    # distances: knn_dists with non-zero items in knn_indices.
    # distances have exactly same number of neighbors for each data point;
    # connectivities have varied number of non-zero items for each data point.
    # connectivities are used as adjacency matrix of the graph for leiden clustering.
    # g = sc._utils.get_igraph_from_adjacency(adata.obsp['connectivities'], directed=True)

    # build knn graph
    n_dims = 20
    n_neighbors = 100
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_dims, use_rep='X_pca')

    df = pd.DataFrame({'cells': adata.obsp['distances'].tocoo().row,
                    'neighbors': adata.obsp['distances'].tocoo().col})
    df = pd.concat([df,
                    pd.DataFrame({'cells': df.cells.unique(),
                                'neighbors': df.cells.unique()})])
    df = df.sort_values(['cells', 'neighbors']).reset_index(drop=True)
    df_knn_indices = df

    # smooth expressions
    assert n_neighbors == df_knn_indices.shape[0]//adata.shape[0]
    indices = df_knn_indices['neighbors'].to_numpy()

    for k, v in markers.items():
        df = pd.DataFrame(index=adata.obs.index)
        for gene in v:
            df['knn'+gene] = adata.raw[:,gene].X[indices.reshape(-1)].reshape(-1,n_neighbors).mean(axis=1)
        adata.obs['knn'+k] = df.mean(axis=1)
    return adata

adata = knn_smooth_expr(adata, markers)

# plot
knn_markers = ['knn' + x for x in markers.keys()]
sc.pl.umap(adata, color=knn_markers, save='_knn_markers.pdf')

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs[0, 0].hist(adata.obs['knnCD3'],50)
axs[0, 0].set_title("knnCD3")
axs[1, 0].hist(adata.obs['knnCD8'],50)
axs[1, 0].set_title("knnCD8")
axs[0, 1].hist(adata.obs['knnCD4'],50)
axs[0, 1].set_title("knnCD4")
axs[1, 1].hist(adata.obs['knnIAC'],50)
axs[1, 1].set_title("knnIAC")
fig.tight_layout()
plt.savefig("./figures/hist_knn_markers.pdf")
plt.clf()

import seaborn as sns
np.random.seed(seed=42)
barcodes_random_sel = np.random.choice(adata.obs.index, size=10000, replace=False)
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
sns.kdeplot(data=adata[barcodes_random_sel].obs, x='knnCD8', y='knnCD4', levels=30, ax=axs[0])
sns.kdeplot(data=adata[barcodes_random_sel].obs, x='knnCD3', y='knnIAC', levels=30, ax=axs[1])
plt.savefig("./figures/kde_knn_markers.pdf")
plt.clf()

# raw expressions
genes = ['CD4', 'CD8A', 'CD8B']
df_raw = pd.DataFrame(adata.raw[:,genes].X.toarray(), columns = genes, index=adata.obs.index)

# classify cell types
cd3_thresh = 0.4
cd8_thresh = 0.9
cd4_thresh = 0.3
iac_thresh = 0.3

adata.obs['cell_type'] = 'Unknown'

idx = (adata.obs['knnCD3'] > cd3_thresh) & ((adata.obs['knnCD8'] > cd8_thresh) | (df_raw['CD8A'] > cd8_thresh) | (df_raw['CD8B'] > cd8_thresh)) & (adata.obs['knnCD4'] < cd4_thresh) & (df_raw['CD4'] < cd4_thresh) 
assert np.alltrue(adata.obs.loc[idx,'cell_type'].to_numpy() == 'Unknown')
adata.obs.loc[idx,'cell_type'] = 'CD8 T'

idx = (adata.obs['knnCD3'] > cd3_thresh) & (adata.obs['knnCD8'] < cd8_thresh) & (df_raw['CD8A'] < cd8_thresh) & (df_raw['CD8B'] < cd8_thresh) & ((adata.obs['knnCD4'] > cd4_thresh) | (df_raw['CD4'] > cd4_thresh))
assert np.alltrue(adata.obs.loc[idx,'cell_type'].to_numpy() == 'Unknown')
adata.obs.loc[idx,'cell_type'] = 'CD4 T'

idx = (adata.obs['knnCD3'] > cd3_thresh) & ((adata.obs['knnCD8'] > cd8_thresh) | (df_raw['CD8A'] > cd8_thresh) | (df_raw['CD8B'] > cd8_thresh)) & ((adata.obs['knnCD4'] > cd4_thresh) | (df_raw['CD4'] > cd4_thresh))
assert np.alltrue(adata.obs.loc[idx,'cell_type'].to_numpy() == 'Unknown')
adata.obs.loc[idx,'cell_type'] = 'DP T'

idx = (adata.obs['knnCD3'] > cd3_thresh) & (adata.obs['knnCD8'] < cd8_thresh) & (adata.obs['knnCD4'] < cd4_thresh) & (df_raw['CD8A'] < cd8_thresh) & (df_raw['CD8B'] < cd8_thresh) & (df_raw['CD4'] < cd4_thresh)
assert np.alltrue(adata.obs.loc[idx,'cell_type'].to_numpy() == 'Unknown')
adata.obs.loc[idx,'cell_type'] = 'DN T'

idx = (adata.obs['knnCD3'] < cd3_thresh) & (adata.obs['knnIAC'] > iac_thresh)
assert np.alltrue(adata.obs.loc[idx,'cell_type'].to_numpy() == 'Unknown')
adata.obs.loc[idx,'cell_type'] = 'IAC'

# finally, as leiden 12 is the IACs-containing cluster, use this cluster info to unambiguously classify IACs and T cells.
adata.obs['cell_type'] = ['Unknown' if ((cluster == '12' and cell_type != 'IAC') or (cluster != '12' and cell_type == 'IAC')) else cell_type for cell_type, cluster in zip(adata.obs.cell_type, adata.obs.leiden)]

# The classified cell types are mostly the same as what we did before, but may not be 100% identical due to slight randomness of embedding/clustering.

# plot
with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": (300)}):
    sc.pl.umap(adata, color='cell_type', size=10, legend_loc='on data', save='_classify_cell_types.pdf')

# write
adata.obs[['cell_type']].to_csv('classify_cell_types.tsv', sep='\t')

