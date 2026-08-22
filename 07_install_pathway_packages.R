# Install Bioconductor if required
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install pathway-analysis packages
BiocManager::install(
  c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
  ),
  ask = FALSE,
  update = FALSE
)

# Test the packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

print("Pathway-analysis packages loaded successfully")