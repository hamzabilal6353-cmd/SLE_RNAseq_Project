
library(SummarizedExperiment)

# Load the expanded dataset
sle_data <- readRDS(
  "data/GSE72509_SLE_recount3_expanded.rds"
)

# Convert sample information into a normal table
metadata_df <- as.data.frame(
  colData(sle_data)
)

# Create Control and SLE labels
metadata_df$condition <- ifelse(
  metadata_df[["sra_attribute.disease_status"]] == "healthy",
  "Control",
  "SLE"
)

# Set Control as the reference group
metadata_df$condition <- factor(
  metadata_df$condition,
  levels = c("Control", "SLE")
)

# Add the condition labels to the RNA-seq dataset
colData(sle_data)$condition <- metadata_df$condition

# Display the number of samples in each group
group_counts <- table(metadata_df$condition)

print("Number of samples in each group:")
print(group_counts)

# Confirm there are 18 controls and 99 SLE patients
stopifnot(
  group_counts["Control"] == 18,
  group_counts["SLE"] == 99
)

# Create a simple metadata table
clean_metadata <- data.frame(
  run_id = as.character(
    metadata_df[["external_id"]]
  ),
  sample_name = as.character(
    metadata_df[["sra.sample_title"]]
  ),
  original_status = as.character(
    metadata_df[["sra_attribute.disease_status"]]
  ),
  condition = as.character(
    metadata_df$condition
  )
)

# Save the clean metadata
write.csv(
  clean_metadata,
  "data/GSE72509_clean_metadata.csv",
  row.names = FALSE
)

# Save the final analysis-ready dataset
saveRDS(
  sle_data,
  "data/GSE72509_SLE_analysis_ready.rds"
)

print("Analysis-ready dataset saved successfully")