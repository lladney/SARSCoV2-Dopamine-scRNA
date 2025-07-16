# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 5: Plot dotplots from GO and KEGG Enrichment

# LIBRARIES
library(clusterProfiler)                                      # clusterProfiler = for GO and KEGG enrichment 
library(ggplot2)                                              # ggplot2 = for dotplot customization of dotplots created with clusterProfiler
library(ggnewscale)                                           # ggnewscale = for multiple colors in ggplot2 plots
library(dplyr)                                                # dplyr = manipulate data (from tidyverse)
library(tibble)                                               # tibble = for tidy data frames (from tidyverse)
library(readr)                                                # readr = for fast file input/output 
library(org.Hs.eg.db)                                         # org.Hs.eg.db = load human gene annotation database for clusterProfiler to map gene IDs and run enrichment analyses

# CREATE RESULTS/FIGURES DIRECTORY 
if (!dir.exists("results/figures")) {                         # Check whether results/figures directory exists
  dir.create("results/figures", recursive = TRUE)             # If directory doesn't exist, create it 
}

# LOADS DEGS TABLE 
degs <- read_csv("results/tables/degs_infected_vs_mock.csv")  # Read CSV file into tibble "degs"

# FILTER FOR SIGNIFICANT GENES
sig_genes <- degs$gene[degs$p_val_adj < 0.05]                 # Filters DEGs for significant ones only, pulling their gene symbols
                                                              # degs$p_val_adj < 0.05 = identifies rows where adjusted p-values < 0.05 (statistically signifcant)
                                                              # degs$gene[...] = extracts gene column of identified statistically significant genes 
                                                              # sig_genes = vector of gene symbols for enrichment
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
go_bp <- enrichGO(gene = genes$ENTREZID,                      # Vector of Entrez IDs to test for enrichment (just filtered these)  
                  OrgDb = org.Hs.eg.db,                       # Use human annotation database for GO term mapping
                  ont = "BP",                                 # Specify ontology: BP = biological process
                  pAdjustMethod = "BH",                       # Adjusts p-values for multiple testing using Benjamin-Hochberg FDR
                  pvalueCutoff = 0.05,                        # Filter out GO terms with raw p-values > 0.05
                  readable = TRUE)                            # Convert Entrez IDs back into gene symbols in outputted result

# KEGG PATHWAY ENRICHMENT FOR BIOLOGICAL PATHWAYS
kegg <- enrichKEGG(gene = genes$ENTREZID,                     # Vector of Entrez IDs to test for enrichment
                   organism = "hsa",                          # Specify organism: hsa = homo sapiens
                   pvalueCutoff = 0.05)                       # Filter to include pathways with raw p-values < 0.05

# SAVE ENRICHMENT RESULTS TO PLOT LATER 
saveRDS(go_bp, file = "results/tables/go_bp.rds")             # Save GO enrichment results to .rds file
saveRDS(kegg,  file = "results/tables/kegg.rds")              # Save KEGG enrichment results to .rds file

# DOTPLOT OF GO BIOLOGICAL PROCESSES
p1 <- dotplot(go_bp, showCategory = 15,                       # Plot top 15 enriched GO terms from go_bp object
              title = "GO Biological Process Enrichment") +   # Dotplot title
  theme_minimal(base_size = 14) +                             # Apply ggplot2 theme with font size 14
  theme(axis.text.y = element_text(size = 10))                # Increase GO term name font size on y-axis

ggsave("results/figures/dotplot_go_bp.png", p1, width = 10, height = 6) # Save dotplot as PNG (10 inches wide x 6 inches tall)

# PLOT KEGG PATHWAYS (IF ENRICHMENT RETURNS RESULTS) 
if (!is.null(kegg) &&                                         # Check if kegg object exists (there are enriched pathways)
    inherits(kegg, "enrichResult") &&                         # Confirms correct object class for plotting
    nrow(as.data.frame(kegg)) > 0) {                          # Confirm at least one enriched pathway
  p2 <- dotplot(kegg, showCategory = 15,                      # Plot top 15 KEGG Pathways from kegg object
                title = "KEGG Pathway Enrichment") +          # Dotplot title
    theme_minimal(base_size = 14) +                           # Apply ggplot2 theme with font size 14
    theme(axis.text.y = element_text(size = 10))              # Increase KEGG Pathway name font size on y-axis
  
  ggsave("results/figures/dotplot_kegg.png", p2, width = 10, height = 6) # Save dotplot as PNG (10 inches wide x 6 inches tall)
} else {                                                      # If there are no enriched pathways...
  message("No KEGG pathways passed the significance threshold.") # Warning message to user to signify no plot generation
}
