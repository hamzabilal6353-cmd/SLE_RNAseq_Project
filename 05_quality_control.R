# Load packages
library(recount3)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

# Create output folders
dir.create("results", showWarnings = FALSE)
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

# Load the analysis-ready dataset
sle_data <- readRDS(
  "data/GSE72509_SLE_analysis_ready.rds"
)

# Transform recount3 coverage counts
assay(sle_data, "counts") <- transform_counts(
  sle_data
)

# Extract transformed counts
count_matrix <- assay(
  sle_data,
  "counts"
)

cat(
  "Original genes:",
  nrow(count_matrix),
  "\n"
)

cat(
  "Samples:",
  ncol(count_matrix),
  "\n"
)

# Keep genes with at least 10 counts in at least 18 samples
keep_genes <- rowSums(
  count_matrix >= 10
) >= 18

filtered_counts <- count_matrix[
  keep_genes,
  ,
  drop = FALSE
]

cat(
  "Genes retained after filtering:",
  nrow(filtered_counts),
  "\n"
)

# Prepare simple sample information
sample_information <- data.frame(
  condition = factor(
    colData(sle_data)$condition,
    levels = c("Control", "SLE")
  ),
  row.names = colnames(sle_data)
)

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_counts),
  colData = sample_information,
  design = ~ condition
)

# Apply variance-stabilising transformation
vst_data <- vst(
  dds,
  blind = TRUE
)

# Calculate PCA
pca_data <- plotPCA(
  vst_data,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar")
)

# Create PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = condition
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  stat_ellipse(
    level = 0.95,
    linewidth = 0.8
  ) +
  xlab(
    paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Control" = "#2471A3",
      "SLE" = "#C0392B"
    )
  ) +
  theme_classic(base_size = 14) +
  ggtitle(
    "PCA of GSE72509: SLE vs Healthy Controls"
  )

# Display the plot
print(pca_plot)

# Save the plot
ggsave(
  filename =
    "results/figures/PCA_SLE_vs_Control.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Save PCA coordinates
write.csv(
  pca_data,
  "results/tables/PCA_sample_coordinates.csv",
  row.names = TRUE
)

# Save objects for the next analysis
saveRDS(
  dds,
  "data/GSE72509_DESeq2_dataset.rds"
)

saveRDS(
  vst_data,
  "data/GSE72509_VST_dataset.rds"
)

cat(
  "PCA plot saved:",
  file.exists(
    "results/figures/PCA_SLE_vs_Control.png"
  ),
  "\n"
)

print("Quality-control step completed")