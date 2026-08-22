# Load required packages
library(recount3)
library(SummarizedExperiment)

# Load the downloaded dataset
sle_data <- readRDS(
  "data/GSE72509_SLE_recount3.rds"
)

# Expand the sample information
sle_data <- expand_sra_attributes(sle_data)

# Extract the metadata table
sample_metadata <- as.data.frame(
  colData(sle_data)
)

# Find useful metadata columns
useful_columns <- grep(
  "title|disease|condition|group|phenotype|attribute",
  colnames(sample_metadata),
  value = TRUE,
  ignore.case = TRUE
)

# Print the column names
cat("Useful metadata columns:\n")
print(useful_columns)

# Print example values from each useful column
for (column_name in useful_columns) {
  cat("\nCOLUMN:", column_name, "\n")
  print(
    head(
      unique(as.character(sample_metadata[[column_name]])),
      20
    )
  )
}

# Save the expanded metadata
write.csv(
  sample_metadata,
  "data/GSE72509_sample_metadata.csv",
  row.names = FALSE
)

# Save the expanded dataset
saveRDS(
  sle_data,
  "data/GSE72509_SLE_recount3_expanded.rds"
)

# Open metadata as a table
View(sample_metadata)