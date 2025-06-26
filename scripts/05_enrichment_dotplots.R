# scripts/05_enrichment_dotplots.R

library(clusterProfiler)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(tibble)
library(readr)
library(org.Hs.eg.db)

# Create results/figures directory if it doesn't exist
if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

# Load DEGs table
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")

# Filter for significant genes
sig_genes <- degs$gene[degs$p_val_adj < 0.05]

# Map to Entrez IDs
genes <- bitr(sig_genes,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
genes <- genes[!is.na(genes$ENTREZID), ]

# Run GO BP enrichment
go_bp <- enrichGO(gene = genes$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# Run KEGG enrichment
kegg <- enrichKEGG(gene = genes$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05)

# Save results
saveRDS(go_bp, file = "results/tables/go_bp.rds")
saveRDS(kegg,  file = "results/tables/kegg.rds")

# Plot GO Biological Processes
p1 <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process Enrichment") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10))
ggsave("results/figures/dotplot_go_bp.png", p1, width = 10, height = 6)

# Plot KEGG Pathways (only if enrichment returned results)
if (!is.null(kegg) && inherits(kegg, "enrichResult") && nrow(as.data.frame(kegg)) > 0) {
  p2 <- dotplot(kegg, showCategory = 15, title = "KEGG Pathway Enrichment") +
    theme_minimal(base_size = 14) +
    theme(axis.text.y = element_text(size = 10))
  ggsave("results/figures/dotplot_kegg.png", p2, width = 10, height = 6)
} else {
  message("No KEGG pathways passed the significance threshold.")
}

