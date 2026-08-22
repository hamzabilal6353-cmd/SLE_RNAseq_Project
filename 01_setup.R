# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install the packages required for the project
BiocManager::install(c(
  "recount3",
  "DESeq2",
  "GEOquery",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "pheatmap"
))

# Install additional plotting packages
install.packages(c(
  "tidyverse",
  "ggrepel",
  "here"
))
# Install Bioconductor package manager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install RNA-seq analysis packages
BiocManager::install(c(
  "recount3",
  "DESeq2",
  "GEOquery",
  "pheatmap"
))

# Install plotting and data-management packages
install.packages(c(
  "tidyverse",
  "ggrepel",
  "here"
))# Install Bioconductor package manager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install RNA-seq analysis packages
BiocManager::install(c(
  "recount3",
  "DESeq2",
  "GEOquery",
  "pheatmap"
))

# Install plotting and data-management packages
install.packages(c(
  "tidyverse",
  "ggrepel",
  "here"
))library(recount3)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(pheatmap)
library(here)

print("All packages loaded successfully")

