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
                                                           # Variable genes stored in: VariableFeatures(seurat)
# SCALE AND RUN PCA
seurat <- ScaleData(seurat)                                # Does 2 things:
                                                           # 1) Centers expression values (subtract mean expression across all cells, so that nean is zero)
                                                           # 2) Scales expression values (divide by standard deviation, so that SD is one)

seurat <- RunPCA(seurat,                                   # Perform Principal Component Analysis (PCA): reduces dimension while preserving variation
                 features = VariableFeatures(seurat))      # Use only those top 2,000 most variable genes identified earlier for PCA (to reduce noise)
                                                           # Note: High-dimensional gene expression data (genes x cells) -> reduced to... principal components (PCs)
                                                           # Note: PC = linear combo of genes (helps explain variance in data), with PC1 being most variant and so on
# VISUALIZE PCA ELBOW PLOT
pdf("results/figures/pca_elbow_plot.pdf")                  # Elbow plot created after this will be saved here
ElbowPlot(seurat)                                          # Generate elbow plot from PCA results (helps determine how many PCs to keep for later steps)
                                                           # Note: PC index (x) by SD (y) to show how much variance explained by each PC
                                                           # Elbow = flat point in curve; before elbow = meaningful variant PCs; after elbow = PCs just noise

dev.off()                                                  # Close the PDF with PCA elbow plot

# RUN UMAP AND CLUSTERING
seurat <- FindNeighbors(seurat, dims = 1:15)               # Use PCs 1-15 (elbow touchdown around 15)
                                                           # Make similarity graph linking cells to their PC neighbors
seurat <- FindClusters(seurat, resolution = 0.5)           # Use similarity fraph from FindNeighbors() to detect clusters; 0.5 = medium resolution
seurat <- RunUMAP(seurat, dims = 1:15)                     # Run Uniform Manifold Approximation and Projection (UMAP) to project the cells into 2D space
                                                           # Note: each point = cell; distance between points = similarity in expression (based on PCA)

# SAVE UMAP PLOT
pdf("results/figures/umap_clusters.pdf")                   # UMAP cluster plot created after this will be saved here
DimPlot(seurat,                                            # Seurat function to plot dimension reductions, one of which is UMAP
        reduction = "umap",                                # Specify UMAP coordinates
        group.by = "seurat_clusters",                      # Color code cells based on cluster assignment
        label = TRUE) +                                    # Add labels to each cluster
  ggtitle("UMAP Clustering of DA Neurons")                 # Add title to the UMAP cluster plot

dev.off()                                                  # Close the PDF with UMAP cluster plot

# SAVE SEURAT OBJECT
saveRDS(seurat, file = "data/processed/combined_clustered.rds") # Save clustered Seurat object in RDS format (avoid repeating clustering in later steps) in processed folder

cat("Normalization, PCA, UMAP, and clustering complete. Output saved to data/processed/combined_clustered.rds\n") # Update the user
