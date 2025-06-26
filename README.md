
# SARSCoV2-Dopamine-scRNA
A modular pipeline for single-cell RNA-seq analysis of SARS-CoV-2-induced gene expression changes in human dopaminergic neurons, featuring Seurat-based clustering, differential expression analysis, and gene ontology enrichment.

## Summary
1. Import of raw expression matrix and quality control using Seurat
2. Normalization, scaling, and clustering of neuronal populations
3. Differential expression analysis (infected vs. mock)
4. Functional enrichment of differentially expressed genes (GO/KEGG)
5. Visualization of enriched pathways

## Project Structure
```
sarscov2-dopamine-scRNA/
├── scripts/                            # R scripts, executed sequentially
│   ├── 01_import_qc.R                  # (1)  Load data, perform initial QC and filtering
│   ├── 02_normalization_clustering.R   # (2)  Normalize, scale, and cluster cells using Seurat
│   ├── 03_differential_expression.R    # (3)  Identify DEGs between infected and mock cells
│   ├── 04_enrichment_analysis.R        # (4)  Perform GO/KEGG enrichment with clusterProfiler
│   ├── 05_enrichment_dotplots.R        # (5)  Generate dotplots for enriched pathways
│   ├── 06_go_enrichment_by_direction.R # (6)  Separate GO terms by up/downregulated DEGs  
│   └── 07_plot_go_dotplots.R           # (7)  Create dotplots highlighting GO categories
│
├── data/
│   └── sample_metadata.csv             # CSV file linking each sample to its condition (infected/mock)
│
├── results/                            # All output files from analysis
│   ├── figures/                        # PDF/PNG plots: UMAPs, PCA, volcano, GO/KEGG dotplots
│   └── tables/                         # CSV tables with DEGs and GO/KEGG enrichment results
│
├── environment.yml                     # Defines Conda environment with R, Seurat, and dependencies
├── .gitignore                          # Excludes unnecessary or large files from Git tracking
└── README.md                           # Project overview, instructions, and structure description

```
## Installation

1. Clone the repository:
```bash
git clone https://github.com/lladney/SARSCoV2-Dopamine-scRNA.git
cd SARSCoV2-Dopamine-scRNA
```  

2. Create a Conda environment:
```bash
conda env create -f environment.yml
conda activate seurat_env
```

## Running the Pipeline

### Step 1:  *IMPORT AND QUALITY CONTROL*
Go to the scripts/ directory and run: 
```r 
Rscript 01_import_qc.R
```
This will: 
- Import the expression matrix
- Filter low-quality cells and genes
- Perform initial UMAP and PCA visualization

### Step 2:  *NORMALIZATION AND CLUSTERING*
Go to the scripts/ directory and run:
```r
Rscript 02_normalization_clustering.R
```
This will: 
- Normalize and scale the data
- Identify variable features
- Perform PCA, clustering, and UMAP dimensionality reduction

### Step 3:  *DIFFERENTIAL EXPRESSION ANALYSIS*
Go to the scripts/ directory and run:
```r
Rscript 03_differential_expression.R
```	
This will: 
- Compare infected vs. mock-treated cells
- Output ```degs_infected_vs_mock.csv``` and volcano plots

### Step 4:  *GO AND KEGG ENRICHMENT ANALYSIS*
Go to the scripts/ directory and run:
```r
Rscript 04_enrichment_analysis.R
```	
This will: 
- Use clusterProfiler to identify enriched GO terms and KEGG pathways
- Output enrichment result tables

### Step 5:  *PLOT ENRICHMENT RESULTS*
Go to the scripts/ directory and run:
```r
Rscript 05_enrichment_dotplots.R
```	
This will: 
- Generate dotplots for enriched pathways

### Step 6:  *GO ENRICHMENT BY DIRECTION*
Go to the scripts/ directory and run:
```r
Rscript 06_go_enrichment_by_direction.R
```	
This will: 
- Separate enrichment results for upregulated vs. downregulated genes

### Step 7:  *VISUALIZE GO TERMS*
Go to the scripts/ directory and run:
```R
Rscript 07_plot_go_dotplots.R
```	
This will: 
- Generate directional dotplots for biological process GO terms

## Notes
* Raw FASTQ files and intermediate large data are excluded via ```.gitignore```
* This pipeline was developed and tested on macOS 10.15 with Conda and R 4.3+
* Expression matrix was derived from GSE248989
