# Required packages for the SLE RNA-seq project

cran_packages <- c(
  "tidyverse",
  "ggrepel",
  "here"
)

bioconductor_packages <- c(
  "recount3",
  "DESeq2",
  "GEOquery",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "pheatmap"
)

# Install BiocManager only when necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing CRAN packages
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# Install missing Bioconductor packages
missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioconductor) > 0) {
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE
  )
}

# Load the main analysis packages
library(recount3)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(here)
library(AnnotationDbi)
library(org.Hs.eg.db)

message("All required packages are installed and loaded successfully.")data_file <- "data/GSE72509_SLE_recount3.rds"

data.frame(
  File = data_file,
  Exists = file.exists(data_file),
  Size_MB = round(file.info(data_file)$size / 1024^2, 2)
)