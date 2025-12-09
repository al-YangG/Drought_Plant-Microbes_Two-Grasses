##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/0_setup_paths.R

source("R/utils_packages.R")
load_project_packages()

# Root is the project directory (assuming you open the .Rproj here)
root <- here::here()

# Define standard subfolders
data_dir    <- file.path(root, "Data")
results_dir <- file.path(root, "results")
fig_dir     <- file.path(root, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)

message("Project root: ", root)
message("Data dir: ", data_dir)
message("Results dir: ", results_dir)
message("Figures dir: ", fig_dir)
