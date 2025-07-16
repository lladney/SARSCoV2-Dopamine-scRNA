# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 4: GO and KEGG enrichment on DEGs (SARS-CoV-2 vs. Mock)

# LIBRARIES
library(tidyverse)                                            # tidyverse = for data wrangling (dplyr), reading/writing files, data frames, exc.       
library(clusterProfiler)                                      # clusterProfiler = for GO and KEGG enrichment 
library(org.Hs.eg.db)                                         # org.Hs.eg.db = load human gene annotation database for clusterProfiler to map gene IDs

# LOAD DEGS TABLE
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")  # Load DEGs table saved in script 03

# FILTER FOR SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES 
sig_genes <- degs %>%                                         # Start tidyverse pipeline with degs data frame
  filter(                                                     # Keep only genes with...
         p_val_adj < 0.05,                                    # 1) adjusted p-value < 0.05 (statistically significant)     
         abs(avg_log2FC) > 0.25) %>%                          # 2) absolute log2 fold change > 0.25 (biologically meaningful change)
  pull(gene)                                                  # Extracts just gene column as character vector -> gives list of gene symbols for enrichment

# CONVERT GENE SYMBOLS TO ENTREZ IDS 
genes <- bitr(                                                # Convert list of gene symbols to Entrez IDs with bitr() function
              sig_genes,                                      # Vector of gene symbols
              fromType = "SYMBOL",                            # Starting with this ID type 
              toType = "ENTREZID",                            # Convert to this ID type (needed for enrichment)
              OrgDb = org.Hs.eg.db)                           # Human annotation database to look up Entrez IDs
                                                              # Result: data frame with two columns, symbol and entrezid
# FILTER OUT UNMAPPED GENES
genes <- genes[!is.na(genes$ENTREZID), ]                      # Remove genes that didn't map to an Entrez ID
                                                              # genes$ENTREZID = mapped Entrez IDs
                                                              # is.na(genes$ENTREZID) = return TRUE if Entrez ID missing
                                                              # ! = invert to keep only mapped genes 
                                                              # genes[... , ] = subset data frame to keep only rows with valid Entrez IDs

# GO ENRICHMENT ANALYSIS FOR BIOLOLGICAL PROCESS
go_bp <- enrichGO(gene         = genes$ENTREZID,              # Vector of Entrez IDs to test for enrichment (just filtered these)         
                  OrgDb        = org.Hs.eg.db,                # Use human annotation database for GO term mapping
                  keyType      = "ENTREZID",                  # Tell function gene IDs are in Entrez format
                  ont          = "BP",                        # Specify ontology: BP = biological process
                                                              # Note: could do other options like MF = molecular function, CG = cellular component
                  pAdjustMethod = "BH",                       # Adjusts p-values for multiple testing using Benjamin-Hochberg FDR
                                                              # Note: controls false discovery rate (FDR) by reducing risk of false positives
                  qvalueCutoff = 0.05,                        # Filter out GO terms with adjusted p-values (q-values) > 0.05, retain only enriched GO terms with FDR < 0.05
                  readable     = TRUE)                        # Convert Entrez IDs back into gene symbols in outputted result

# KEGG PATHWAY ENRICHMENT FOR BIOLOGICAL PATHWAYS
kegg <- enrichKEGG(gene         = genes$ENTREZID,             # Vector of Entrez IDs to test for enrichment
                   organism     = "hsa",                      # Specify organism: hsa = homo sapiens
                   pAdjustMethod = "BH",                      # Adjusts p-values for multiple testing using Benjamin-Hochberg FDR
                   qvalueCutoff = 0.05)                       # Filter to include pathways with adjusted p-values (FDR) < 0.05

# SAVE RESULTS
write_csv(as_tibble(go_bp), "results/tables/go_enrichment.csv")  # Save enriched biological process terms in DEGs to csv 
write_csv(as_tibble(kegg), "results/tables/kegg_enrichment.csv") # Save enriched biological pathways in DEGs to csv

# PRINT TOP RESULTS
print(head(go_bp, 10))                                        # Print first 10 rows of GO enrichment results (top BP terms) 
print(head(kegg, 10))                                         # Print first 10 rows of KEGG pathway results (top pathways) 

# SAVE ENRICHMENT RESULTS TO PLOT LATER 
saveRDS(go_bp, file = "results/tables/go_bp.rds")             # Save GO enrichment results to .rds file
saveRDS(kegg, file = "results/tables/kegg.rds")               # Save KEGG enrichment results to .rds file

cat("Enrichment analysis complete. Results saved to results/tables/\n") # Update the user
