# Load required packages
library(recount3)
library(SummarizedExperiment)

# Load the saved SLE dataset
sle_data <- readRDS(
  "data/GSE72509_SLE_recount3.rds"
)

# Show the number of genes and samples
print("Number of genes and samples:")
print(dim(sle_data))

# Show available count-data types
print("Available assay names:")
print(assayNames(sle_data))

# Extract sample information
sample_metadata <- as.data.frame(
  colData(sle_data)
)

# Show metadata dimensions
print("Metadata dimensions:")
print(dim(sample_metadata))

# Find columns that may contain patient/control information
metadata_columns <- grep(
  "title|sample|condition|disease|group|character|attribute",
  colnames(sample_metadata),
  value = TRUE,
  ignore.case = TRUE
)

print("Potential patient/control columns:")
print(metadata_columns)

# Save all sample information as a CSV file
write.csv(
  sample_metadata,
  file = "data/GSE72509_sample_metadata.csv",
  row.names = FALSE
)

# Open the sample information in RStudio
View(sample_metadata)

print("Dataset inspection completed")