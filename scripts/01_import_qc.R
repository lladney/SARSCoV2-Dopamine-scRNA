# 01_import_qc.R
# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons (GSE248989)
# Step 1: Import, QC, merge Seurat objects

# ---- Libraries ----
library(Seurat)
library(tidyverse)
library(patchwork)

# ---- Paths ----
data_dirs <- list(
  "infected" = "data/raw/infected/",
  "mock"     = "data/raw/mock/"
)

# Metadata file
metadata_path <- "data/metadata/sample_metadata.csv"

# ---- Load sample metadata ----
sample_metadata <- read_csv(metadata_path)

# ---- Read and create Seurat objects ----
seurat_list <- map2(
  .x = data_dirs,
  .y = names(data_dirs),
  .f = function(path, sample_id) {
    counts <- Read10X(data.dir = path)
    obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
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
