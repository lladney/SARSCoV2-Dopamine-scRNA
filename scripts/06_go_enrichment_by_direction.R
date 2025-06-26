# scripts/06_go_enrichment_by_direction.R

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# Load DEGs
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")

# Filter for significant DEGs
sig_degs <- degs %>% filter(p_val_adj < 0.05)

# Separate upregulated and downregulated
up_genes <- sig_degs %>% filter(avg_log2FC > 0) %>% pull(gene)
down_genes <- sig_degs %>% filter(avg_log2FC < 0) %>% pull(gene)

# Convert to ENTREZ IDs
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Identify unmapped symbols
unmapped_up <- setdiff(up_genes, up_entrez$SYMBOL)
unmapped_down <- setdiff(down_genes, down_entrez$SYMBOL)

# Combine into a data frame
unmapped_df <- data.frame(
  gene_symbol = c(unmapped_up, unmapped_down),
  regulation = c(rep("upregulated", length(unmapped_up)),
                 rep("downregulated", length(unmapped_down)))
)

# Save to CSV
write.csv(unmapped_df, "results/tables/unmapped_genes.csv", row.names = FALSE)

# Print to terminal
message(glue::glue(
  "{length(unmapped_df$gene_symbol)} unmapped gene symbols saved to: results/tables/unmapped_genes.csv\n",
  "   - Upregulated: {length(unmapped_up)}\n",
  "   - Downregulated: {length(unmapped_down)}"
))

# Filter out NAs
up_entrez <- up_entrez[!is.na(up_entrez$ENTREZID), ]
down_entrez <- down_entrez[!is.na(down_entrez$ENTREZID), ]

# GO enrichment for upregulated genes
go_up <- enrichGO(gene = up_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# GO enrichment for downregulated genes
go_down <- enrichGO(gene = down_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

# Save results
saveRDS(go_up, file = "results/tables/go_up.rds")
saveRDS(go_down, file = "results/tables/go_down.rds")

write.csv(as.data.frame(go_up), file = "results/tables/go_up.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), file = "results/tables/go_down.csv", row.names = FALSE)

