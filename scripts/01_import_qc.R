# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 1: Import, QC, merge Seurat objects

# LIBRARIES
library(Seurat)                                            # Seurat = for scRNAseq analysis
library(tidyverse)                                         # Tidyverse = for data manipulation, visualization, analysis
library(patchwork)                                         # Patchwork = for combining multiple ggplot2 plots

# PATHS
data_dirs <- list(                                         # Create list for directory paths
  "infected" = "data/raw/infected/",                       # Path to COVID-infected neuronal data
  "mock"     = "data/raw/mock/"                            # Path to mock-infected neuronal data
)

metadata_path <- "data/metadata/sample_metadata.csv"       # Define file path string and store in variable

sample_metadata <- read_csv(metadata_path)                 # Use readr package (tidyverse) to read in metadata CSV and store in data frame
                                                           # sample_metadata object -> tibble (tidy data frame)
# READ AND CREATE SEURAT OBJECTS
seurat_list <- map2(                                       # Use purrr package (tidyverse) to iterate over parallel vectors .x and .y and...
                                                           # ...apply function .f to each pair of elements
                                                           # Results of each iteration returned as list
  .x = data_dirs,                                          # Loop through .x (paths to data)
  .y = names(data_dirs),                                   # Loop through .y (names: infected and mock)
  .f = function(path, sample_id) {                         # Function reads over gene expression matrix from each path, creates Seurat object,...
                                                           # ...annotates it with sample ID and condition from metadata, and returns Seurat object
                                                           # Executes twice (once for infected, then mock)
    
    counts <- Read10X(data.dir = path)                     # Read in gene expression count data from 10x Genomics Cell Ranger pipeline
                                                           # Looks for: 
                                                           # 1) matrix.mtx.gz (gene-by-cell sparse matrix)
                                                           # 2) features.tsv.gz (list of features: genes, IDs)
                                                           # 3) barcodes.tsv.gz (cell barcode identifiers)
                                                           # Returns: sparse matrix R object with rows (genes), columns (cells), entries = UMI counts
   
    obj <- CreateSeuratObject(                             # Raw gene expression matrix -> Seurat object for scRNAseq analysis
                              counts = counts,             # Gene (row) by cell (column) sparse dgCMatrix with UMI couts in each cell, read from Read10x()
                              project = sample_id,         # Assign project name to Seurat object (sample_id = infected or mock)
                              min.cells = 3,               # Filter out genes detected in <3 cells (remove rare genes)
                              min.features = 200)          # Filter out cells detected with <200 genes (standard to filter out low-quality cells... 
                                                           # ...i.e., empty droplets, dead cells, lysed cells) 
                                                           # Could adjust down since working with neurons which have lower RNA content
    
    obj$sample_id <- sample_id                             # Add new column sample_id (infected or mock) to metadata of Seurat object
    obj$condition <- sample_metadata %>%                   # Add new column condition to metadata of Seurat object
    filter(sample_id == !!sample_id) %>%                   # Keep only rows where sample_id matches current sample ID being processed
    pull(condition)                                        # Use pull() from dplyr (tidyverse) to extract condition from matched row into vector
    return(obj)                                            # Return Seurat object defined in map2()
                                                           # obj contains: 
                                                           # 1) UMI count matrix (genes x cells)
                                                           # 2) Initial QC metrics (nFeature_RNA, nCount_RNA)
                                                           # 3) Two metadata fields: a) sample_id: mock or infected, and b) condition: control or infected
  }
)

# MERGE DATASETS
combined <- merge(                                         # Merge two Seurat objects (infected, mock) into combined object (combined)
                  seurat_list$infected,                    # First Seurat object in merge
                  y = seurat_list$mock,                    # Second Seurat object to merge into the first
                  add.cell.ids = names(seurat_list),       # Add prefix to cell barcode (indicates original sample) so barcodes are unique
                  project = "SARS-CoV-2_DA")               # Set project name for new combined object

# CALCULATE MITOCHONDRIAL GENE PERCENTAGE
combined[["percent.mt"]] <-                                # Store percentage in Seurat object's metadata under percent.mt
PercentageFeatureSet(combined, pattern = "^MT-")           # Calculate proportion of transcripts (UMI counts) from (mitochondrial) genes starting with MT-
                                                           # Note: mitochondrial genes = highly expressed when cell is under high stress
# VISUALIZE QC METRICS  
VlnPlot(combined,                                          # Create violin plots for metadata features in Seurat object 
                                                           # Note: violin plots show distribution of values across all cells for a feature
        features = c(                                      # Plot the following (QC metrics)...
  "nFeature_RNA",                                          # 1) Number of genes (features) detected per cell
  "nCount_RNA",                                            # 2) Total number of UMI counts per cell
  "percent.mt"),                                           # 3) Percent of reads mapping to mitochondrial genes
        ncol = 3)                                          # Arrangement: 1 row, 3 columns
                                                           # Result: 3 violin plots
                                                           # 1) nFeature_RNA: filters out low-complexity cells or doublets (when high)
                                                           # 2) nCount_RNA: detects outliers (when high)
                                                           # 3) percent.mt: flags stressed or lysed cells (when high)
                                                           # Note: Low everything = empty droplets

# FILTER CELLS
combined <- subset(combined,                               # subset() returns new Seurat object with filtered cells only
                                                           # Cells must meet these criteria to stay:
                   subset = nFeature_RNA > 300 &           # 1) Cells express fewer than 300 genes? Out. Likely empty droplets or dead cells
                            nFeature_RNA < 5000 &          # 2) Cells express more than 5000 genes? Out. Likely doublets containing multiple cells in one droplet
                            percent.mt < 10)               # 3) More than 10% mitochondrial content? Out. Likely stressed/lysed/dying cells
                                                           # Now just filtered, cleaned cells remain
# SAVE PREPROCESSED DATA
saveRDS(combined, file = "data/processed/combined_qc.rds") # Save filtered Seurat object in RDS format (avoid preprocessing in later steps) in processed folder
                                                           # RDS file will include: count matrix, metadata, QC metrics (percent.mt, sample_id, condition)
                                                           # Load in future: combined <- readRDS("data/processed/combined_qc.rds")

cat("Import and QC complete. Output saved to data/processed/combined_qc.rds\n") # Update the user
