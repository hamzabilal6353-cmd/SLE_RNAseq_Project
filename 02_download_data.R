# Load recount3
library(recount3)

# Allow extra time for downloading
options(timeout = 1200)

# Create the data folder
dir.create("data", showWarnings = FALSE)

# Get available human RNA-seq projects
human_projects <- available_projects()

# Find the SLE project
project_info <- subset(
  human_projects,
  project == "SRP062966" &
    project_type == "data_sources"
)

# Show the project information
print(project_info)

# Confirm that the project was found
if (nrow(project_info) == 0) {
  stop("SRP062966 was not found in recount3.")
}

# Download the gene-count dataset
sle_data <- create_rse(project_info)

# Save the downloaded dataset
saveRDS(
  sle_data,
  file = "data/GSE72509_SLE_recount3.rds"
)

# Show the number of genes and samples
print(dim(sle_data))

print("SLE RNA-seq data downloaded and saved successfully")