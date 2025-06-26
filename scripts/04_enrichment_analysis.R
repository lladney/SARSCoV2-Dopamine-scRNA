# 04_enrichment_analysis.R
# GO and KEGG enrichment on DEGs (SARS-CoV-2 vs. Mock)
# Load required libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load DEGs table
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")

# Filter for significantly differentially expressed genes
sig_genes <- degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
genes <- bitr(sig_genes,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)

# Filter out unmapped genes
genes <- genes[!is.na(genes$ENTREZID), ]

# GO enrichment analysis (Biological Process)
go_bp <- enrichGO(gene         = genes$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable     = TRUE)

# KEGG pathway enrichment (optional)
kegg <- enrichKEGG(gene         = genes$ENTREZID,
                   organism     = "hsa",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

# Save results
write_csv(as_tibble(go_bp), "results/tables/go_enrichment.csv")
write_csv(as_tibble(kegg), "results/tables/kegg_enrichment.csv")

# Print top results
print(head(go_bp, 10))
print(head(kegg, 10))

# Save enrichment results for downstream plotting
saveRDS(go_bp, file = "results/tables/go_bp.rds")
saveRDS(kegg, file = "results/tables/kegg.rds")

cat("Enrichment analysis complete. Results saved to results/tables/\n")