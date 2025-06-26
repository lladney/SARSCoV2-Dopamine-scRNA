
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
│
├── .gitignore                          # Excludes unnecessary or large files from Git tracking
│
└── README.md                           # Project overview, instructions, and structure description

```
## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/plasmodium-rnaseq-pipeline.git
cd plasmodium-rnaseq-pipeline
```  

2. Create a Conda environment:
```bash
conda create -n plasmodium_env python=3.10
conda activate plasmodium_env
```

3. Install Python dependencies:
```
pip install -r requirements.txt
```

4. Install the required command-line tools via [bioconda](https://bioconda.github.io/):
   - `fastqc`: Quality control for raw reads
   - `cutadapt`: Adapter trimming
   - `salmon`: Transcript-level quantification
   - `multiqc`: Aggregated QC reporting

   You can install them all at once:
   ```bash
   conda install -c bioconda fastqc cutadapt salmon multiqc

## Running the Pipeline

### Step 1:  *PREPROCESSING*
Go to the preprocessing/ directory and run: 
```python 
Preprocessing_Pipeline_1.py
```
This will perform: 
- GEO/SRA Metadata Fetching
- FASTQ Download
- Adapter Trimming with Cutadapt
- Quality Control with FastQC
- MultiQC Summary

### Step 2:  *QUANTIFICATION*
Go to the quantification/ directory and run:
```python
Alignment_Quantification_2.py
```
This will: 
- Build a Salmon Index
- Quantify transcript expression
- Generate gene-level count matrix and metadata CSV

### Step 3:  *DIFFERENTIAL EXPRESSION ANALYSIS*
Go to the dgeanalysis/ directory and run:
```Rscript
DGE_Analysis_DESeq2.R
```	
This will produce: 
- deseq2_results.csv: Differential expression output
- volcano_plot.png: Volcano plot of DEGs
- heatmap_top20_DEGs.png: Heatmap of top 20 DEGs
- PCA_plot.png: Principal component analysis of samples

## Notes
* Raw data files and large intermediate results are .gitignored
* This pipeline was developed and tested on macOS 10.15 with Conda and R 4.3+
