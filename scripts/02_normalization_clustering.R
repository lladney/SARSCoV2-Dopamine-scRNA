# 02_normalization_clustering.R
# Normalize, scale, reduce dimensions, and cluster SARS-CoV-2 vs. mock cells

library(Seurat)
library(ggplot2)

# Load QC'd object
seurat <- readRDS("data/processed/combined_qc.rds")

# Normalize and identify variable features
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

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