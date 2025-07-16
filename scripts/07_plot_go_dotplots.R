# SARS-CoV-2 scRNA-seq in hPSC-derived DA neurons
# Step 7: Create Dotplots for Up- vs. Down-regulated Genes

# LIBRARIES
library(clusterProfiler)                                        # clusterProfiler = for interpreting and plotting GO results 
library(ggplot2)                                                # ggplot2 = create plots (from tidyverse)

# LOAD GO ENRICHMENT RESULTS
go_up <- readRDS("results/tables/go_up.rds")                    # Read in upregulated genes objects from .rds files created in script 06
go_down <- readRDS("results/tables/go_down.rds")                # Read in downregulated genes objects from .rds files created in script 06

# CREATE DOTPLOT FOR UPREGULATED GENES  
p_up <- dotplot(go_up, showCategory = 15,                       # Plot top 15 enriched GO terms from go_up object 
                title = "GO Enrichment: Upregulated Genes") +   # Dotplot title
  theme_minimal(base_size = 14) +                               # Apply ggplot2 theme with font size 14
  theme(axis.text.y = element_text(size = 10))                  # Increase GO term name font size on y-axis

# CREATE DOTPLOT FOR DOWNREGULATED GENES 
p_down <- dotplot(go_down, showCategory = 15,                   # Plot top 15 enriched GO terms from go_down object 
                  title = "GO Enrichment: Downregulated Genes") + # Dotplot title
  theme_minimal(base_size = 14) +                               # Apply ggplot2 theme with font size 14
  theme(axis.text.y = element_text(size = 10))                  # Increase GO term name font size on y-axis

# SAVE THE PLOTS 
ggsave("results/figures/dotplot_go_upregulated.png", p_up, width = 10, height = 6)      # Save upregulated dotplot as PNG (10 inches wide x 6 inches tall)
ggsave("results/figures/dotplot_go_downregulated.png", p_down, width = 10, height = 6)  # Save downregulated dotplot as PNG (10 inches wide x 6 inches tall)
