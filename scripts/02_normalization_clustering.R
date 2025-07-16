# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons (GSE248989)
# Step 2: Normalize, scale, reduce dimensions, and cluster COVID vs mock cells

# LIBRARIES
library(Seurat)                                            # Seurat = for scRNAseq analysis
library(ggplot2)                                           # ggplot2 = create plots (from tidyverse)

# LOAD PREPROCESSED SEURAT OBJECT
seurat <- readRDS("data/processed/combined_qc.rds")        # Load filtered Seurat object saved in script 01

# NORMALIZE AND IDENTIFY VARIABLE FEATURES
seurat <- NormalizeData(seurat)                            # Normalize gene expression for each cell
                                                           # Note: Important because cells have high range of UMI counts so normalization makes cell comparable
seurat <- FindVariableFeatures(seurat,                     # Identifies most variable genes across all cells
                               selection.method = "vst",   # vst = Variance Stabilizing Transformation, which models gene expression mean/variance 
                                                           # Note: preferred for high-throughput scRNAseq data 
                               nfeatures = 2000)           # Select top 2,000 most variable genes
                                                           # Note: likely show biological variation, not just noise

# Scale and run PCA
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

# Visualize PCA elbow plot
pdf("results/figures/pca_elbow_plot.pdf")
ElbowPlot(seurat)
dev.off()

# Run UMAP and clustering
seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.5)
seurat <- RunUMAP(seurat, dims = 1:15)

# Save UMAP plot
pdf("results/figures/umap_clusters.pdf")
DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP Clustering of DA Neurons")
dev.off()

# Save Seurat object
saveRDS(seurat, file = "data/processed/combined_clustered.rds")

cat("Normalization, PCA, UMAP, and clustering complete. Output saved to data/processed/combined_clustered.rds\n")
