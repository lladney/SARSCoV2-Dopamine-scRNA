# scripts/07_plot_go_dotplots.R

library(clusterProfiler)
library(ggplot2)

# Load GO enrichment results
go_up <- readRDS("results/tables/go_up.rds")
go_down <- readRDS("results/tables/go_down.rds")

# Create dotplot for upregulated genes

p_up <- dotplot(go_up, showCategory = 15, title = "GO Enrichment: Upregulated Genes") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10))

# Create dotplot for downregulated genes
p_down <- dotplot(go_down, showCategory = 15, title = "GO Enrichment: Downregulated Genes") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10))

# Save the plots
ggsave("results/figures/dotplot_go_upregulated.png", p_up, width = 10, height = 6)
ggsave("results/figures/dotplot_go_downregulated.png", p_down, width = 10, height = 6)


