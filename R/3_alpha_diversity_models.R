##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/3_alpha_diversity_models.R
##
## Tasks:
##  1) Fit LMMs for microbial richness:
##       Richness ~ Moisture * Climate * Species + (1 | Mesocosm)
##       Richness ~ GWC      * Climate * Species + (1 | Mesocosm)
##       Richness ~ VWC      * Climate * Species + (1 | Mesocosm)
##  2) Export Type-III ANOVA tables for Moisture, GWC, VWC
##  3) Compute variance explained (%) for Moisture models
##     using variance_explained_lmm() from utils_stats.R

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_stats.R")

load_project_packages()

root        <- here::here()
data_dir    <- file.path(root, "Data", "3_Alpha diversity")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Import alpha-diversity tables (per compartment)
# -------------------------------------------------------------------

div_BR <- read.csv(file.path(data_dir, "Alpha Diversity_BR.csv"),
                   check.names = FALSE)
div_BS <- read.csv(file.path(data_dir, "Alpha Diversity_BS.csv"),
                   check.names = FALSE)
div_FL <- read.csv(file.path(data_dir, "Alpha Diversity_FL.csv"),
                   check.names = FALSE)
div_FR <- read.csv(file.path(data_dir, "Alpha Diversity_FR.csv"),
                   check.names = FALSE)
div_FS <- read.csv(file.path(data_dir, "Alpha Diversity_FS.csv"),
                   check.names = FALSE)

div_list <- list(
  BR = div_BR,  # Bacteria - Root endosphere
  BS = div_BS,  # Bacteria - Rhizosphere
  FL = div_FL,  # Fungi    - Leaf endosphere
  FR = div_FR,  # Fungi    - Root endosphere
  FS = div_FS   # Fungi    - Rhizosphere
)

# Ensure consistent factor structure
div_list <- lapply(div_list, function(dat) {
  dat$Climate  <- factor(dat$Climate)
  dat$Species  <- factor(dat$Species)
  dat$Mesocosm <- factor(dat$Mesocosm)
  dat
})

# -------------------------------------------------------------------
# 2) Helper: run Richness LMMs for a given continuous predictor
# -------------------------------------------------------------------

run_richness_models <- function(div_list, cont_var) {
  # Fit LMMs: Richness ~ cont_var * Climate * Species + (1 | Mesocosm)
  models <- lapply(div_list, function(dat) {
    fit_lmm_3way(dat          = dat, response     = "Richness",
                 cont_var     = cont_var, random_factor = "Mesocosm")
  })
  
  # Type-III ANOVA for each compartment
  anova_list <- lapply(models, anova_type3)
  
  # Combine ANOVA tables and add FDR-adjusted p-values
  anova_tbl_list <- mapply(
    FUN = function(a, comp) {
      df <- as.data.frame(a)
      df$Term <- rownames(df)
      rownames(df) <- NULL
      
      # identify p-value column (Pr(>F) or Pr(>Chisq))
      p_col <- grep("^Pr\\(>", names(df), value = TRUE)[1]
      
      df$Compartment <- comp
      df$p_adj_FDR   <- p.adjust(df[[p_col]], method = "fdr")
      
      df[, c("Compartment", "Term",
             setdiff(names(df), c("Compartment", "Term")))]
    },
    a    = anova_list,
    comp = names(anova_list),
    SIMPLIFY = FALSE
  )
  
  anova_tbl <- do.call(rbind, anova_tbl_list)
  
  list(models = models, anova = anova_list, anova_tbl = anova_tbl)
}

# -------------------------------------------------------------------
# 3) Richness ~ Moisture × Climate × Species
# -------------------------------------------------------------------

res_moisture <- run_richness_models(div_list, cont_var = "Moisture")

# 3.1 Save ANOVA table (Moisture)
write.csv(res_moisture$anova_tbl,
          file.path(results_dir, "ANOVA_Moisture.csv"), ## Table 2
          row.names = FALSE)

# 3.2 Variance explained (%) for Moisture models
var_exp_moisture_list <- mapply(
  FUN = function(a, mod, comp) {
    df_var <- variance_explained_lmm(a, mod)  # from utils_stats.R
    df_var$Compartment <- comp
    df_var
  },
  a    = res_moisture$anova,
  mod  = res_moisture$models,
  comp = names(res_moisture$models),
  SIMPLIFY = FALSE
)

var_exp_moisture_tbl <- do.call(rbind, var_exp_moisture_list)

write.csv(var_exp_moisture_tbl,
          file.path(results_dir, "VarExp_Richness.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 4) Richness ~ GWC × Climate × Species
# -------------------------------------------------------------------

res_gwc <- run_richness_models(div_list, cont_var = "GWC")

write.csv(res_gwc$anova_tbl,
          file.path(results_dir, "ANOVA_GWC.csv"), ## Table S6
          row.names = FALSE)

# -------------------------------------------------------------------
# 5) Richness ~ VWC × Climate × Species
# -------------------------------------------------------------------

res_vwc <- run_richness_models(div_list, cont_var = "VWC")

write.csv(res_vwc$anova_tbl,
          file.path(results_dir, "ANOVA_VWC.csv"), row.names = FALSE)

