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

# signatures
# Cytotoxicity: 
sig_cyto = ["CX3CR1", "FCGR3A", "FGFBP2", "GNLY", "GZMA", "GZMB", "GZMH", "KLRG1", "NKG7", "PRF1"]
# Dysfunction: 
sig_dysfun = ["CD274", "PDCD1", "LAG3", "HAVCR2", "CTLA4", "TIGIT", "CD96", "KLRC1", "CD160", "TIMD4", "CD27", "TNFRSF9", "ICOS", "CD28", "TNFRSF18", "TMIGD2", "TNFSF14", "CD226", "TNFRSF4", "TNFRSF8"]
# Memory (CD8 T): 
sig_cd8tmem = ["LEF1",  "TCF7", "SERINC5", "IL7R", "CCR7", "TNFRSF25", "S1PR1", "PASK", "FLT3LG", "CAMK4", "SORL1", "DHRS3", "TMEM63A", "MGAT4A", "LTB", "RCAN3", "ABLIM1", "PLAC8", "DGKA", "TC2N", "SELL", "KLRB1", "NOL4L", "TESPA1", "MCUB", "GIMAP5", "OXNAD1", "FAM102A", "SATB1", "NOSIP", "RIPOR2", "ICAM2", "ATM", "SCML4", "CD5", "PIK3IP1", "FOXP1", "EPB41", "CD28", "GOLGA8A", "GIMAP7", "CHMP7", "NELL2", "DENND2D", "GOLGA8B", "GIMAP2", "TMEM123", "GPR183", "TTC39C", "TMEM131L", "RAPGEF6", "AAK1"]

gs = {'sig_cyto': sig_cyto, 
      'sig_dysfun': sig_dysfun,
      'sig_cd8tmem': sig_cd8tmem}

# calculate sig scores
for sig_name, sig_genes in gs.items():
    sc.tl.score_genes(adata, sig_genes, score_name=sig_name)
adata.obs[gs.keys()].to_csv('sig_cyto_dysfun_cd8mem.tsv', sep='\t')

sc.pl.umap(adata, color=gs.keys(), save='_sig_cyto_dysfun_cd8mem.pdf')

