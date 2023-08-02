#!/usr/bin/env Rscript

rm(list = ls())

library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ggpubr)


args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Need passing arguments to the R script from command lines")

sample_path <- args[1]

data_dir <- paste0(sample_path, 
                   "/filtered_feature_bc_matrix/")


pdf("doublets.pdf", height = 4, width = 5)


# load data and create seurat object -----
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
n_cells_recovered <- ncol(data)

# keep all data so that each one has a score/classification; lysed cells may contribute to doublets
seu_obj <- CreateSeuratObject(counts = data)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)

# Perform dimensionality reduction by PCA and UMAP embedding 
seu_obj <- RunPCA(seu_obj, npcs = 100, verbose = FALSE)
n_dims <- 100
seu_obj <- RunUMAP(seu_obj, dims = 1:n_dims, verbose = FALSE)
seu_obj <- FindNeighbors(seu_obj, dims = 1:n_dims, verbose = FALSE)
seu_obj <- FindClusters(seu_obj, resolution = 0.6, verbose = FALSE)


# doublet scores -------
sweep.res.list_seu_obj <- paramSweep_v3(seu_obj, PCs = 1:n_dims, sct = FALSE) # no SCTransform
sweep.stats_seu_obj <- summarizeSweep(sweep.res.list_seu_obj, GT = FALSE)
bcmvn_seu_obj <- find.pK(sweep.stats_seu_obj)

bcmvn_seu_obj %>%
  as_tibble() %>%
  write_tsv("bcmvn.tsv")

id_best <- which.max(bcmvn_seu_obj$BCmetric)
pK_best <- levels(bcmvn_seu_obj$pK)[bcmvn_seu_obj$pK[id_best]] %>% as.double()

# emprical multiplet rate (derived from 10X multiplet rate table)
get_rate <- function(n_cells_recovered) {
  if (n_cells_recovered > 1000) {
    return((0.8 + (n_cells_recovered - 1000)*0.8/1000)*0.01)
  } else {
    return((0.4 + (n_cells_recovered - 500)*0.4/500)*0.01)
  }
}

## Homotypic Doublet Proportion Estimate 
doublet_rate <- round(get_rate(n_cells_recovered), digits = 3)
annotations <- seu_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_obj@meta.data$ClusteringResults
nExp_poi <- round(doublet_rate*nrow(seu_obj@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:n_dims, pN = 0.25, pK = pK_best, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# rerun using adjusted rate
pANN_name <- paste0("pANN_0.25_", pK_best, "_", nExp_poi)
seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:n_dims, pN = 0.25, pK = pK_best, nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = FALSE)

# scatter of score
FeaturePlot(seu_obj, pANN_name, reduction = "pca")
FeaturePlot(seu_obj, pANN_name)

# score vs count
plot1 <- ggscatter(seu_obj@meta.data, "nCount_RNA", pANN_name)
plot2 <- ggscatter(seu_obj@meta.data, "nFeature_RNA", pANN_name)
plot1 + plot2

# doublet class based on adjusted rate
df_detectable_rate <- round(nExp_poi.adj/nrow(seu_obj@meta.data), digits = 3)
cls_adj_name <- paste0("DF.classifications_0.25_", pK_best, "_", nExp_poi.adj)
DimPlot(seu_obj, reduction = "pca", group.by = cls_adj_name) + 
  ggtitle(paste0("detected rate: ", df_detectable_rate))
DimPlot(seu_obj, group.by = cls_adj_name)  + 
  ggtitle(paste0("detected rate: ", df_detectable_rate))

# doublet class vs count
plot1 <- ggboxplot(seu_obj@meta.data, cls_adj_name, "nCount_RNA")
plot2 <- ggboxplot(seu_obj@meta.data, cls_adj_name, "nFeature_RNA")
plot1 + plot2

# write df scores
tb <- seu_obj@meta.data[, (colnames(seu_obj@meta.data) %in% c("nCount_RNA", "nFeature_RNA", "seurat_clusters")) |
                          grepl("^pANN", colnames(seu_obj@meta.data)) |
                          grepl("^DF", colnames(seu_obj@meta.data))] %>%
  as.tibble()
tibble(barcode = rownames(seu_obj@meta.data)) %>%
  bind_cols(tb) %>%
  write_tsv("doublets.tsv")


dev.off()

