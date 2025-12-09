##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/1_plant_performance.R
## Plant performance: biomass, proportion, SLA, LDMC
## Reproduces plant models.

source("R/utils_packages.R")
source("R/utils_stats.R")

load_project_packages()

root       <- here::here()
data_dir   <- file.path(root, "Data", "1_Plant performance")
results_dir <- file.path(root, "results")
fig_dir     <- file.path(root, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)

# --- Load data ---------------------------------------------------------------

data_biomass <- read.csv(file.path(data_dir, "Species biomass.csv"), header = TRUE)
data_traits  <- read.csv(file.path(data_dir, "SLA+LDMC.csv"), header = TRUE)

# Documentation of key columns is in the original script, not repeated here.

# --- Models: plant biomass & proportion -------------------------------------

# Linear models for species biomass
lm_moisture <- lm(Species_biomass ~ Moisture * Climate * Species,
                  data = data_biomass)
lm_gwc      <- lm(Species_biomass ~ GWC * Climate * Species,
                  data = data_biomass)

res_lm_moisture <- car::Anova(lm_moisture, type = "II")
res_lm_gwc      <- car::Anova(lm_gwc,      type = "II")

# Beta regression for species proportion (0–100%, convert to 0–1)
data_biomass <- data_biomass |>
  dplyr::mutate(Proportion_adj = Proportion * 0.01)

beta_moisture <- betareg::betareg(Proportion_adj ~ Moisture * Climate * Species,
                                  data = data_biomass)

beta_gwc <- betareg::betareg(Proportion_adj ~ GWC * Climate * Species,
                             data = data_biomass)

res_beta_moisture <- car::Anova(beta_moisture, type = "III")
res_beta_gwc      <- car::Anova(beta_gwc,      type = "III")

# --- Traits: SLA & LDMC (LMMs) ----------------------------------------------

data_traits <- data_traits |>
  dplyr::mutate(Mesocosm = factor(Mesocosm))

sla_moisture  <- lme4::lmer(SLA  ~ Moisture * Climate * Species + (1 | Mesocosm),
                            data = data_traits)
ldmc_moisture <- lme4::lmer(LDMC ~ Moisture * Climate * Species + (1 | Mesocosm),
                            data = data_traits)

sla_gwc  <- lme4::lmer(SLA  ~ GWC * Climate * Species + (1 | Mesocosm),
                       data = data_traits)
ldmc_gwc <- lme4::lmer(LDMC ~ GWC * Climate * Species + (1 | Mesocosm),
                       data = data_traits)

res_sla_moisture  <- car::Anova(sla_moisture,  type = "II")
res_ldmc_moisture <- car::Anova(ldmc_moisture, type = "II")
res_sla_gwc       <- car::Anova(sla_gwc,       type = "II")
res_ldmc_gwc      <- car::Anova(ldmc_gwc,      type = "II")

# --- Combine Moisture models into one table ---------------------------------

as_df_std <- function(tab, response, model_type) {
  df <- as.data.frame(tab)
  df <- std_names_anova(df)
  df$Response  <- response
  df$ModelType <- model_type
  
  common_cols <- c("Response", "ModelType", "Df", "Statistic", "p_value")
  df[, intersect(common_cols, names(df)), drop = FALSE]
}

df_lm_m   <- as_df_std(res_lm_moisture,   "Species biomass",   "LM")
df_beta_m <- as_df_std(res_beta_moisture, "Biomass proportion","Beta")
df_sla_m  <- as_df_std(res_sla_moisture,  "SLA",               "LMM")
df_ldmc_m <- as_df_std(res_ldmc_moisture, "LDMC",              "LMM")

table_moisture <- dplyr::bind_rows(df_lm_m, df_beta_m, df_sla_m, df_ldmc_m)

write.csv(table_moisture,
          file.path(results_dir, "Table_Moisture_AllResults.csv"),
          row.names = TRUE) ## Table 1

# --- Combine GWC models -----------------------------------------------------

df_lm_g   <- as_df_std(res_lm_gwc,   "Species biomass",   "LM")
df_beta_g <- as_df_std(res_beta_gwc, "Biomass proportion","Beta")
df_sla_g  <- as_df_std(res_sla_gwc,  "SLA",               "LMM")
df_ldmc_g <- as_df_std(res_ldmc_gwc, "LDMC",              "LMM")

table_gwc <- dplyr::bind_rows(df_lm_g, df_beta_g, df_sla_g, df_ldmc_g)

write.csv(table_gwc,
          file.path(results_dir, "Table_GWC_AllResults.csv"),
          row.names = TRUE) ## Table S5

# --- Variance explained (%) for Moisture-based models (plant traits) --------

# 1) Species biomass (LM)
df_var_lm <- {
  df <- as.data.frame(res_lm_moisture)
  df$Term <- rownames(df)
  df <- df[df$Term != "Residuals", ]
  SS_total <- sum(df[["Sum Sq"]], na.rm = TRUE)
  R2 <- summary(lm_moisture)$r.squared
  df$Variance_explained_percent <- 100 * (df[["Sum Sq"]] / SS_total) * R2
  df$Variance_explained_percent <- round(df$Variance_explained_percent, 4)
  df$Response <- "Species biomass"
  df[, c("Response", "Term", "Variance_explained_percent")]
}

# 2) SLA (LMM)
df_var_sla <- variance_explained_lmm(res_sla_moisture, sla_moisture)
# 3) LDMC (LMM)
df_var_ldmc <- variance_explained_lmm(res_ldmc_moisture, ldmc_moisture)

table_var_exp_plant <- dplyr::bind_rows(df_var_lm, df_var_sla, df_var_ldmc) |>
  dplyr::arrange(Response, Term)

write.csv(table_var_exp_plant,
          file.path(results_dir, "Variance_explained_plant.csv"),
          row.names = FALSE)
