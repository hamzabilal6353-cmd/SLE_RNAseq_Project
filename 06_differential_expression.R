# Load packages
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)

# Load the DESeq2 dataset
dds <- readRDS(
  "data/GSE72509_DESeq2_dataset.rds"
)

# Run differential-expression analysis
dds <- DESeq(dds)

# Compare SLE against healthy controls
de_results <- results(
  dds,
  contrast = c(
    "condition",
    "SLE",
    "Control"
  ),
  alpha = 0.05
)

# Convert results into a normal table
results_table <- as.data.frame(
  de_results
)

# Add Ensembl identifiers
results_table$ensembl_id <- rownames(
  results_table
)

# Load original dataset to retrieve gene symbols
sle_data <- readRDS(
  "data/GSE72509_SLE_analysis_ready.rds"
)

gene_information <- as.data.frame(
  rowData(sle_data)
)

# Add gene symbols when available
if ("gene_name" %in% colnames(gene_information)) {
  results_table$gene_symbol <-
    gene_information[
      match(
        results_table$ensembl_id,
        rownames(gene_information)
      ),
      "gene_name"
    ]
} else {
  results_table$gene_symbol <- NA_character_
}

# Classify genes
results_table$classification <- "Not significant"

results_table$classification[
  !is.na(results_table$padj) &
    results_table$padj < 0.05 &
    results_table$log2FoldChange >= 1
] <- "Upregulated"

results_table$classification[
  !is.na(results_table$padj) &
    results_table$padj < 0.05 &
    results_table$log2FoldChange <= -1
] <- "Downregulated"

# Sort by adjusted p-value
results_table <- results_table[
  order(results_table$padj),
]

# Create significant-results tables
significant_results <- subset(
  results_table,
  classification != "Not significant"
)

upregulated_results <- subset(
  results_table,
  classification == "Upregulated"
)

downregulated_results <- subset(
  results_table,
  classification == "Downregulated"
)

# Display result counts
cat(
  "Genes tested:",
  nrow(results_table),
  "\n"
)

cat(
  "Significant genes:",
  nrow(significant_results),
  "\n"
)

cat(
  "Upregulated genes:",
  nrow(upregulated_results),
  "\n"
)

cat(
  "Downregulated genes:",
  nrow(downregulated_results),
  "\n"
)

# Save result tables
write.csv(
  results_table,
  "results/tables/DESeq2_all_genes.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "results/tables/DESeq2_significant_genes.csv",
  row.names = FALSE
)

write.csv(
  upregulated_results,
  "results/tables/DESeq2_upregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  downregulated_results,
  "results/tables/DESeq2_downregulated_genes.csv",
  row.names = FALSE
)

# Prepare data for volcano plot
volcano_data <- subset(
  results_table,
  !is.na(padj)
)

volcano_data$minus_log10_padj <- -log10(
  pmax(volcano_data$padj, 1e-300)
)

# Select the ten most statistically significant genes
label_data <- head(
  volcano_data[
    order(volcano_data$padj),
  ],
  10
)

label_data$plot_label <- ifelse(
  is.na(label_data$gene_symbol) |
    label_data$gene_symbol == "",
  label_data$ensembl_id,
  label_data$gene_symbol
)

# Create volcano plot
volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    colour = classification
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = plot_label),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = c(
      "Downregulated" = "#2471A3",
      "Not significant" = "grey70",
      "Upregulated" = "#C0392B"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Differential Expression: SLE vs Control",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    colour = "Classification"
  )

print(volcano_plot)

# Save volcano plot
ggsave(
  "results/figures/Volcano_SLE_vs_Control.png",
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# Save fitted DESeq2 dataset
saveRDS(
  dds,
  "data/GSE72509_DESeq2_fitted.rds"
)

print("Differential-expression analysis completed")