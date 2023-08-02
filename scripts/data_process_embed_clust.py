import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import re

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=300, facecolor='white')
# sc.settings.figdir = './figures/' # default

results_file = 'ana.h5ad'  # the file that will store the analysis results

# download counts from GEO and put them into ../data with each directory corresponding to a sample
# load sample_ids
df = pd.read_csv("../data/samples.tsv", sep="\t")
df['path'] = ['../data/'+ x for x in df.sample_id]
df_sample = df

# read and concat data
adatas = []
for p in df_sample.path:
    print(p)
    adata = sc.read_10x_mtx(
        p,  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        )
    adatas.append(adata)

adata = adatas[0].concatenate(*adatas[1:],
    batch_key="sample_id",
    uns_merge="unique",
    batch_categories=df_sample.sample_id)

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

# load doublets info predicted by DoubletFinder
df = pd.read_csv('../data/doublets_allsamples.tsv', sep='\t')
df.barcode = df.barcode + '-' + df.sample_id
df = df.drop(['sample_id'], axis=1)

# rm doublets
barcodes_sel = df[df.DoubletFinder_class != "Doublet"].barcode.str.replace("^(.*)_(.*)$", "\\2-\\1", regex=True).to_list()
adata = adata[adata.obs_names.isin(barcodes_sel)]

# Preprocessing
sc.pl.highest_expr_genes(adata, n_top=20, save='.pdf')

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# A violin plot of some of the computed quality measures:
# the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_genes_counts_mt.pdf')

# Remove cells that have too many mitochondrial genes expressed or too many total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_counts_mt.pdf')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_counts_genes.pdf')
sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', save='_genes_mt.pdf')

adata = adata[adata.obs.n_genes_by_counts < 7000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]

# Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# hvg
# It is necessary for hvg to include variable genes across batches so that the embedding has batch variations, which are to be corrected using Harmony
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=5000) #, batch_key='sample_id')
sc.pl.highly_variable_genes(adata, save='.pdf')

# tcr genes
df = pd.read_csv("../data/hgnc_gene_group_tcr.txt", sep="\t")
df_ex = df

tcr_ig_genes = df_ex['Approved symbol'].to_list()

# intersect of hvg and tcr genes
len(set(adata.var_names[adata.var.highly_variable]) & set(tcr_ig_genes))

# mask tcr, pick top hvg, mask the rest of hvg
adata.var.loc[adata.var_names.isin(tcr_ig_genes), 'highly_variable'] = False

n_hvg = 3000
hvg_sorted = adata.var[adata.var.highly_variable].highly_variable_rank.sort_values()[:n_hvg].index
adata.var.loc[~adata.var_names.isin(hvg_sorted), 'highly_variable'] = False
adata.var.loc[~adata.var_names.isin(hvg_sorted), 'highly_variable_rank'] = np.nan

# Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use.
adata.raw = adata

# do the filtering
adata = adata[:, adata.var.highly_variable]

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. 
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

# Load cell cycle genes defined in Tirosh et al, 2015. It is a list of 97 genes, represented by their gene symbol. The list here is for humans, in case of alternate organism, a list of ortologues should be compiled. There are major differences in the way Scanpy and Seurat manage data, in particular we need to filter out cell cycle genes that are not present in our dataset to avoid errors.
cell_cycle_genes = [x.strip() for x in open('../data/regev_lab_cell_cycle_genes.txt')]

# Here we define two lists, genes associated to the S phase and genes associated to the G2M phase
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.raw.var_names]

# We here perform cell cycle scoring. The function is actually a wrapper to sc.tl.score_gene_list, which is launched twice, to score separately S and G2M phases. Both sc.tl.score_gene_list and sc.tl.score_cell_cycle_genes are a port from Seurat and are supposed to work in a very similar way. To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list. Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Now we can regress out both S score and G2M score.
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

# Principal component analysis
n_pcs = 100
sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, save='.pdf')
sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, save='_logscale.pdf')

# pca load
sc.pl.pca_loadings(adata, components = range(1,11), save='.pdf')

# plot cumsum of var ratio
df = pd.DataFrame.from_dict(
    {"PC": list(range(1, (n_pcs+1))),
    "variance_ratio": adata.uns["pca"]["variance_ratio"]})
df["variance_ratio_cumsum"] = df.variance_ratio.cumsum()
df.plot(x="PC", y="variance_ratio_cumsum", style='.-')
plt.savefig("./figures/pca_variance_ratio_cumsum.pdf")

# determine dim
n_dims = 20
print(n_dims, 'PCs explain a total variance ratio of', adata.uns["pca"]["variance_ratio"][:n_dims].sum())

# run harmony
sce.pp.harmony_integrate(adata, 'sample_id', max_iter_harmony=100)

# Computing the neighborhood graph
sc.pp.neighbors(adata, n_pcs=n_dims, use_rep='X_pca_harmony')

# Embedding the neighborhood graph
sc.tl.umap(adata)

# Clustering the neighborhood graph
sc.tl.leiden(adata, resolution=0.7)

with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": (300)}):
    sc.pl.umap(adata, color='leiden', size=10, legend_loc='on data', save='_leiden.pdf')

with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": (300)}):
    sc.pl.umap(adata, color='sample_id', size=10, save='_batch.pdf')

with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": (300)}):
    sc.pl.umap(adata, color='phase', size=10, save='_cell_cycle.pdf')

# The results are mostly the same as what we did before, but may not be 100% identical due to slight randomness of embedding/clustering (e.g., package version, order of samples loaded in the counts matrix).
# As the IACs-containing cluster also include T cells, to resolve the heterogeneity, we further sub-clustered the IACs-containing cluster and got the fine-grained cluster of IACs.

# Write
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

