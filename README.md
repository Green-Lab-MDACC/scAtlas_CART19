# Single cell atlas of CD19 CAR T-cells

We have described a resource of single cell RNA-sequencing (scRNA-seq) data from the infusion products of 59 relapsed or refractory large B-cell lymphoma (rrLBCL) patients treated with axicabtagene ciloleucel (Axi-cel; 417,167 single cells), that we have made publicly available to facilitate ongoing discovery efforts that will improve patient outcomes.

Here we provide the scripts and data for processing the scRNA-seq counts data to reproduce the resource.

## Code Description
- `data_process_embed_clust.py` processes the scRNA-seq counts data of all samples to generate an AnnData object of all cells.
- `classify_cell_types.py` classifies all cells into CD8+, CD4+ T cells, and other cell types.
- `sig_cyto_dysfun_cd8mem.py` calculates scores of cytotoxicity, dysfunction, and memory signatures.
- `doublets.R` uses DoubletFinder to predict doublets/singlets for cells of individual samples.

