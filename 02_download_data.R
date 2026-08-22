# Load recount3
library(recount3)

# Allow extra time for downloading
options(timeout = 1200)

# Create a folder for the data
dir.create("data", showWarnings = FALSE)

# Get the list of available human RNA-seq projects
human_projects <- available_projects()

# Find the SLE study: SRP062966 / GSE72509
project_info <- subset(
  human_projects,
  project == "SRP062966" &
    project_type == "data_sources"
)

# Display the matching project
print(project_info)

# Stop if the study was not found
if (nrow(project_info) == 0) {
  stop("SRP062966 was not found in recount3.")
}

# Download and create the gene-count dataset
sle_data <- create_rse(project_info)

# Save it so it does not need to be downloaded again
saveRDS(
  sle_data,
  file = "data/GSE72509_SLE_recount3.rds"
)

# Display the number of genes and samples
print(dim(sle_data))

print("SLE RNA-seq data downloaded and saved successfully")
