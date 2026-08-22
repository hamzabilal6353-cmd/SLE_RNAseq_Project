# Load packages
library(DESeq2)
library(SummarizedExperiment)
library(pheatmap)

# Create output folders if necessary
dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Load variance-stabilised expression data
vst_data <- readRDS(
  "data/GSE72509_VST_dataset.rds"
)

# Read differential-expression results
de_results <- read.csv(
  "results/tables/DESeq2_all_genes.csv",
  stringsAsFactors = FALSE
)

# Extract expression matrix
vst_matrix <- SummarizedExperiment::assay(
  vst_data
)

# Retain significant genes present in the expression matrix
significant_results <- subset(
  de_results,
  !is.na(padj) &
    classification != "Not significant" &
    ensembl_id %in% rownames(vst_matrix)
)

# Sort genes by adjusted p-value
significant_results <- significant_results[
  order(significant_results$padj),
]

# Remove duplicated identifiers
significant_results <- significant_results[
  !duplicated(significant_results$ensembl_id),
]

# Stop if no genes were found
if (nrow(significant_results) == 0) {
  stop("No significant genes were available for the heatmap.")
}

# Select the top 30 genes
number_of_genes <- min(
  30,
  nrow(significant_results)
)

top_genes <- head(
  significant_results,
  number_of_genes
)

# Extract expression values for selected genes
heatmap_matrix <- vst_matrix[
  top_genes$ensembl_id,
  ,
  drop = FALSE
]

# Create readable gene labels
gene_labels <- ifelse(
  is.na(top_genes$gene_symbol) |
    top_genes$gene_symbol == "",
  top_genes$ensembl_id,
  top_genes$gene_symbol
)

# Ensure that every row label is unique
gene_labels <- make.unique(
  gene_labels
)

rownames(heatmap_matrix) <- gene_labels

# Convert expression values into row Z-scores
heatmap_matrix_z <- t(
  scale(
    t(heatmap_matrix)
  )
)

# Replace any non-finite values with zero
heatmap_matrix_z[
  !is.finite(heatmap_matrix_z)
] <- 0

# Prepare sample-condition annotation
sample_condition <- as.character(
  SummarizedExperiment::colData(
    vst_data
  )$condition
)

sample_annotation <- data.frame(
  Condition = factor(
    sample_condition,
    levels = c(
      "Control",
      "SLE"
    )
  )
)

rownames(sample_annotation) <- colnames(
  heatmap_matrix_z
)

# Set sample annotation colours
annotation_colours <- list(
  Condition = c(
    "Control" = "#2471A3",
    "SLE" = "#C0392B"
  )
)

# Set heatmap colours
heatmap_colours <- colorRampPalette(
  c(
    "#2166AC",
    "#F7F7F7",
    "#B2182B"
  )
)(100)

# Set the displayed Z-score range
heatmap_breaks <- seq(
  -3,
  3,
  length.out = 101
)

# Save the heatmap as PNG
pheatmap(
  heatmap_matrix_z,
  scale = "none",
  color = heatmap_colours,
  breaks = heatmap_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_col = sample_annotation,
  annotation_colors = annotation_colours,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  main = "Top 30 Differentially Expressed Genes: SLE vs Control",
  filename = "results/figures/Top30_DEG_heatmap.png",
  width = 14,
  height = 10
)

# Save a PDF version
pheatmap(
  heatmap_matrix_z,
  scale = "none",
  color = heatmap_colours,
  breaks = heatmap_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_col = sample_annotation,
  annotation_colors = annotation_colours,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  main = "Top 30 Differentially Expressed Genes: SLE vs Control",
  filename = "results/figures/Top30_DEG_heatmap.pdf",
  width = 14,
  height = 10
)

# Save the selected gene information
write.csv(
  top_genes[
    ,
    c(
      "ensembl_id",
      "gene_symbol",
      "log2FoldChange",
      "padj",
      "classification"
    )
  ],
  "results/tables/Top30_DEG_heatmap_genes.csv",
  row.names = FALSE
)

# Display completion information
cat(
  "Genes included in heatmap:",
  nrow(heatmap_matrix_z),
  "\n"
)

cat(
  "Samples included in heatmap:",
  ncol(heatmap_matrix_z),
  "\n"
)

cat(
  "Heatmap saved:",
  file.exists(
    "results/figures/Top30_DEG_heatmap.png"
  ),
  "\n"
)

print("Top-gene heatmap completed successfully")