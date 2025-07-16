# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 3: Identify DEGs between SARS-CoV-2 and mock DA neurons

# LIBRARIES
library(Seurat)                                            # Seurat = for scRNAseq analysis
library(ggplot2)                                           # ggplot2 = create plots (from tidyverse)
library(dplyr)                                             # dplyr = manipulate data (from tidyverse)
library(tibble)                                            # tibble = for tidy data frames (from tidyverse)

# LOAD CLUSTERED SEURAT OBJECT
seurat <- readRDS("data/processed/combined_clustered.rds") # Load clustered Seurat object saved in script 02

# RECOVER CONDITION INFO IF MISSING
if (is.null(seurat$condition)) {                           # Check if condition column exists in Seurat object's metadata
                                                           # If TRUE, column didn't transfer over from previous step
  seurat$condition <- ifelse(                              # Reconstruct condition from cell names (barcodes)
    grepl("^infected_",                                    # Check which barcodes start with "infected" -> assume COVID-infected samples
          colnames(seurat)), 
    "SARS-CoV-2",                                          # Label cells with "SARS-CoV-2" if infected
    "Mock")                                                # If that cell doesn't have "infected" prefix in barcode, labeled as "Mock"
}

# SET IDENTITY
Idents(seurat) <- seurat$condition                         # Treat each cell's condition (infected or mock) as main grouping
seurat <- JoinLayers(seurat)                               # Ensure Seurat object runs on full dataset (2 merged 10x datasets, infected and mock)

# Run differential expression
de_markers <- FindMarkers(seurat, ident.1 = "SARS-CoV-2", ident.2 = "Mock", 
                          logfc.threshold = 0.25, min.pct = 0.1)

# Add gene column and sort
de_markers <- de_markers %>% 
  rownames_to_column("gene") %>%
  arrange(p_val_adj)

# Save to CSV
write.csv(de_markers, "results/tables/degs_infected_vs_mock.csv", row.names = FALSE)

# Basic volcano plot
pdf("results/figures/volcano_degs.pdf")
ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "DEGs: SARS-CoV-2 vs. Mock",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")
dev.off()

cat("DE analysis complete. Output saved to results/tables/degs_infected_vs_mock.csv\n")
