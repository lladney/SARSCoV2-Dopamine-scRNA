# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 6: Perform GO Enrichment for Up- vs. Down-regulated Genes

# LIBRARIES
library(clusterProfiler)                                         # clusterProfiler = running enrichment analyses
library(org.Hs.eg.db)                                            # org.Hs.eg.db = human genome annotation database for gene ID mapping
library(dplyr)                                                   # dplyr = manipulate data (from tidyverse)
library(readr)                                                   # readr = fast CSV file input/output

# LOAD DEGS
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")     # Load DEGs table saved in script 03

# FILTER FOR SIGNIFICANT DEGS
sig_degs <- degs %>% filter(p_val_adj < 0.05)                    # Filter out genes with adjusted p-value > 0.05

# SEPARATE UP- AND DOWN-REGULATED GENES
up_genes <- sig_degs %>%                                         # From significant DEGs (sig_degs)...
  filter(avg_log2FC > 0) %>%                                     # Select genes with positive log2 fold change -> upregulated in COVID condition
  pull(gene)                                                     # Extract gene symbol column

down_genes <- sig_degs %>%                                       # From significant DEGs (sig_degs)...
  filter(avg_log2FC < 0) %>%                                     # Select genes with negative log2 fold change -> downregulated in COVID condition
  pull(gene)                                                     # Extract gene symbol column

# CONVERT TO ENTREZ IDS
up_entrez <- bitr(up_genes,                                      # Use bitr() from clusterProfiler to convert upregulated gene symbols to Entrez IDs
                  fromType = "SYMBOL",                           # Starting with gene symbol ID type
                  toType = "ENTREZID",                           # Convert to Entrez ID
                  OrgDb = org.Hs.eg.db)                          # Human annotation database to look up Entrez IDs

down_entrez <- bitr(down_genes,                                  # Use bitr() from clusterProfiler to convert downregulated gene symbols to Entrez IDs
                    fromType = "SYMBOL",                         # Starting with gene symbol ID type
                    toType = "ENTREZID",                         # Convert to Entrez ID
                    OrgDb = org.Hs.eg.db)                        # Human annotation database to look up Entrez IDs

# IDENTIFY UNMAPPED GENE SYMBOLS
unmapped_up <- setdiff(up_genes, up_entrez$SYMBOL)               # Find gene symbols in up_genes not in mapped result (up_entrez$SYMBOL)
unmapped_down <- setdiff(down_genes, down_entrez$SYMBOL)         # Find gene symbols in down_genes not in mapped result (down_entrez$SYMBOL)

# COMBINE INTO DATA FRAME
unmapped_df <- data.frame(                                       # Organize unmapped gene symbols into tidy data frame
  gene_symbol = c(unmapped_up, unmapped_down),                   # Concatenates all unmapped gene symbols into one column
  regulation = c(rep("upregulated", length(unmapped_up)),        # Create second column indicating whether gene is up- or down-regulated
                 rep("downregulated", length(unmapped_down)))    # Second column will contain label "upregulated" or "downregulated"
)

# SAVE TO CSV
write.csv(unmapped_df, "results/tables/unmapped_genes.csv",      # Save unmapped_df data frame to CSV file
          row.names = FALSE)                                     # Stop R from writing row column number into CSV

# PRINT # GENES FAILED TO MAP TO TERMINAL
message(glue::glue(
  "{length(unmapped_df$gene_symbol)} unmapped gene symbols saved to: results/tables/unmapped_genes.csv\n",
  "   - Upregulated: {length(unmapped_up)}\n",
  "   - Downregulated: {length(unmapped_down)}"
))

# FILTER OUT GENE SYMBOLS WITHOUT ENTREZ ID
up_entrez <- up_entrez[!is.na(up_entrez$ENTREZID), ]             # Keep only rows of upregulated gene symbols where Entrez ID is not missing
down_entrez <- down_entrez[!is.na(down_entrez$ENTREZID), ]       # Keep only rows of downregulated gene symbols where Entrez ID is not missing

# GO ENRICHMENT FOR UPREGULATED GENES
go_up <- enrichGO(                                               # Identify upregulated GO terms significantly enriched
                  gene = up_entrez$ENTREZID,                     # Vector of Entrez IDs for upregulated genes
                  OrgDb = org.Hs.eg.db,                          # Human gene annotation database
                  keyType = "ENTREZID",                          # Specify ID type (Entrez ID)
                  ont = "BP",                                    # GO ontology for Biological Process
                  pAdjustMethod = "BH",                          # Benjamini–Hochberg FDR correction
                  pvalueCutoff = 0.05,                           # Filter to include BPs with raw p-values < 0.05
                  qvalueCutoff = 0.2)                            # Filter to include BPs with adjusted p-values (FDR) < 0.2
  
# GO ENRICHMENT FOR DOWNREGULATED GENES
go_down <- enrichGO(                                             # Identify downregulated GO terms significantly enriched
                    gene = down_entrez$ENTREZID,                 # Vector of Entrez IDs for downregulated genes
                    OrgDb = org.Hs.eg.db,                        # Human gene annotation database
                    keyType = "ENTREZID",                        # Specify ID type (Entrez ID)
                    ont = "BP",                                  # GO ontology for Biological Process
                    pAdjustMethod = "BH",                        # Benjamini–Hochberg FDR correction
                    pvalueCutoff = 0.05,                         # Filter to include BPs with raw p-values < 0.05
                    qvalueCutoff = 0.2)                          # Filter to include BPs with adjusted p-values (FDR) < 0.2

# SAVE RESULTS TO .RDS FILES 
saveRDS(go_up, file = "results/tables/go_up.rds")                # Save upregulated GO enrichment results to .rds file
saveRDS(go_down, file = "results/tables/go_down.rds")            # Save downregulated GO enrichment results to .rds file

write.csv(as.data.frame(go_up), file = "results/tables/go_up.csv", row.names = FALSE)       # Save upregulated enriched biological process terms in DEGs to csv 
write.csv(as.data.frame(go_down), file = "results/tables/go_down.csv", row.names = FALSE)   # Save downregulated enriched biological process terms in DEGs to csv
