# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons (GSE248989)
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
    obj$sample_id <- sample_id
    obj$condition <- sample_metadata %>% filter(sample_id == !!sample_id) %>% pull(condition)
    return(obj)
  }
)

# ---- Merge datasets ----
combined <- merge(seurat_list$infected, y = seurat_list$mock, add.cell.ids = names(seurat_list), project = "SARS-CoV-2_DA")

# ---- Calculate mitochondrial gene percentage ----
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# ---- Visualize QC metrics ----
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ---- Filter cells ----
combined <- subset(combined,
                   subset = nFeature_RNA > 300 & 
                            nFeature_RNA < 5000 & 
                            percent.mt < 10)

# ---- Save output ----
saveRDS(combined, file = "data/processed/combined_qc.rds")

cat("Import and QC complete. Output saved to data/processed/combined_qc.rds\n")
