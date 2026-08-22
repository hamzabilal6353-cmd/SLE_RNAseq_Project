# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Create output folders
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Read differential-expression results
de_results <- read.csv(
  "results/tables/DESeq2_all_genes.csv",
  stringsAsFactors = FALSE
)

# Remove version numbers from Ensembl identifiers
de_results$ensembl_clean <- sub(
  "\\..*$",
  "",
  de_results$ensembl_id
)

# Convert Ensembl identifiers to Entrez identifiers
gene_mapping <- bitr(
  unique(de_results$ensembl_clean),
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# Prepare the background of all tested genes
background_genes <- unique(
  gene_mapping$ENTREZID
)

# Identify upregulated Ensembl genes
upregulated_ensembl <- unique(
  de_results$ensembl_clean[
    de_results$classification == "Upregulated"
  ]
)

# Identify downregulated Ensembl genes
downregulated_ensembl <- unique(
  de_results$ensembl_clean[
    de_results$classification == "Downregulated"
  ]
)

# Convert upregulated genes to Entrez identifiers
upregulated_entrez <- unique(
  gene_mapping$ENTREZID[
    gene_mapping$ENSEMBL %in% upregulated_ensembl
  ]
)

# Convert downregulated genes to Entrez identifiers
downregulated_entrez <- unique(
  gene_mapping$ENTREZID[
    gene_mapping$ENSEMBL %in% downregulated_ensembl
  ]
)

# GO Biological Process analysis for upregulated genes
go_upregulated <- enrichGO(
  gene = upregulated_entrez,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# GO Biological Process analysis for downregulated genes
go_downregulated <- enrichGO(
  gene = downregulated_entrez,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Convert results to normal tables
go_up_table <- as.data.frame(
  go_upregulated
)

go_down_table <- as.data.frame(
  go_downregulated
)

# Save pathway tables
write.csv(
  go_up_table,
  "results/tables/GO_BP_upregulated.csv",
  row.names = FALSE
)

write.csv(
  go_down_table,
  "results/tables/GO_BP_downregulated.csv",
  row.names = FALSE
)

# Create and save upregulated pathway plot
if (nrow(go_up_table) > 0) {
  
  go_up_plot <- enrichplot::dotplot(
    go_upregulated,
    showCategory = 15
  ) +
    ggtitle(
      "GO Biological Processes: Upregulated Genes"
    )
  
  print(go_up_plot)
  
  ggsave(
    "results/figures/GO_BP_upregulated.png",
    plot = go_up_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Create and save downregulated pathway plot
if (nrow(go_down_table) > 0) {
  
  go_down_plot <- enrichplot::dotplot(
    go_downregulated,
    showCategory = 15
  ) +
    ggtitle(
      "GO Biological Processes: Downregulated Genes"
    )
  
  print(go_down_plot)
  
  ggsave(
    "results/figures/GO_BP_downregulated.png",
    plot = go_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# Save pathway-analysis objects
saveRDS(
  go_upregulated,
  "results/tables/GO_BP_upregulated_object.rds"
)

saveRDS(
  go_downregulated,
  "results/tables/GO_BP_downregulated_object.rds"
)

# Display summary
cat(
  "Mapped background genes:",
  length(background_genes),
  "\n"
)

cat(
  "Mapped upregulated genes:",
  length(upregulated_entrez),
  "\n"
)

cat(
  "Mapped downregulated genes:",
  length(downregulated_entrez),
  "\n"
)

cat(
  "Significant upregulated pathways:",
  nrow(go_up_table),
  "\n"
)

cat(
  "Significant downregulated pathways:",
  nrow(go_down_table),
  "\n"
)

print("GO pathway analysis completed successfully")