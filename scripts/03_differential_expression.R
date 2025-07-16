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

# RUN DIFFERENTIAL EXPRESSION TO IDENTIFY UP- VS DOWN-REGULATED GENES  
de_markers <- FindMarkers(seurat,                          # Compare two cell groups, which are...
                          ident.1 = "SARS-CoV-2",          # 1) The test group (infected cells)
                          ident.2 = "Mock",                # 2) The reference group (mock cells)
                          logfc.threshold = 0.25,          # Filter for genes with ~20% fold change in either direction
                          min.pct = 0.1)                   # Filter for genes expressed in >10% of cells in either group

# ADD GENE COLUMN AND SORT DE RESULTS
de_markers <- de_markers %>%                               # Start pipe to apply transformations to de_markers data frame
  rownames_to_column("gene") %>%                           # Converts: gene names as row names -> row names into new column "gene"
  arrange(p_val_adj)                                       # Sorts: most to least statistically significant (top to bottom)

# SAVE TO CSV
write.csv(de_markers,                                      # Write de_markers data frame to CSV file
          "results/tables/degs_infected_vs_mock.csv",      # File path and name
          row.names = FALSE)                               # Tell R not to write extra index column for row numbers (already moved gene names to "gene" column

# BASIC VOLCANO PLOT
pdf("results/figures/volcano_degs.pdf")                    # Volcano plot created after this will be saved here
ggplot(de_markers,                                         # Initialize volcano plot with de_markers data frame
       aes(x = avg_log2FC,                                 # Map to x-axis: log2 fold change between infected and mock
                                                           # Note: genes on... right = upregulated, left = downregulated
           y = -log10(p_val_adj))) +                       # Map to y-axis: inverse log scale of adjusted p-value
                                                           # Note: the higher the gene, the more significant
           geom_point(alpha = 0.5) +                       # Add transparent point for each gene (based on log fold change and p-value)
                                          
           geom_hline(yintercept = -log10(0.05),           # Draw red horizontal dashed line at -log10(0.05) to mark the adjusted p-value significance threshold
                      linetype = "dashed",                 # Note: genes above line? signficant; below? not significant
                      color = "red") +
           geom_vline(xintercept = c(-1, 1),               # Draw two vertical blue dashed lines at -1 and 1 to mark 2-fold change in either direction
                      linetype = "dashed",                 # Note: highlights genes of strong biological effect sizes
                      color = "blue") +
           theme_minimal() +                               # Apply clean background
           
          labs(                                            # Add labels
           title = "DEGs: SARS-CoV-2 vs. Mock",           
           x = "Log2 Fold Change",                        
           y = "-Log10 Adjusted P-Value")

dev.off()                                                  # Close the PDF with volcano plot

cat("DE analysis complete. Output saved to results/tables/degs_infected_vs_mock.csv\n") # Update the user
