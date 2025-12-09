##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## R/utils_packages.R

load_project_packages <- function() {
  pkgs <- c(
    # data handling
    "dplyr", "tidyr", "purrr", "tibble", "readr", "data.table",
    # modelling
    "lme4", "car", "betareg", "MASS", "MuMIn", "emmeans", "multcomp",
    # community ecology
    "phyloseq", "vegan", "MicEco", "permute",
    # plotting
    "ggplot2", "ggpubr", "patchwork", "ggplotify", "ggrepel",
    # misc
    "here"
  )
  
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing package: ", p)
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
  
  # Type-III sums-of-squares by default
  options(contrasts = c("contr.sum", "contr.poly"))
}


