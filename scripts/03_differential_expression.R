# 03_differential_expression.R
# Identify DEGs between SARS-CoV-2 and mock DA neurons

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)

# Load clustered object
seurat <- readRDS("data/processed/combined_clustered.rds")

# Recover condition info if missing
if (is.null(seurat$condition)) {
  seurat$condition <- ifelse(grepl("^infected_", colnames(seurat)), "SARS-CoV-2", "Mock")
}

# Set identity
Idents(seurat) <- seurat$condition
seurat <- JoinLayers(seurat)

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
