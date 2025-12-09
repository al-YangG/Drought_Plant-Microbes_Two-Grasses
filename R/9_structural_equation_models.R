##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: October 2025 ##-##
##
## scripts/9_structural_equation_models.R
##
## Tasks:
##   1) Import SEM-ready dataset (traits, microbial metrics, biomass)
##   2) Compute species-specific correlation matrices (Festuca, Lolium)
##   3) Visualize significant correlations (upper-triangle heatmap)
##   4) Fit species-specific SEMs (piecewiseSEM)
##   5) Extract standardized coefficients and R² for comparison
##   6) Export correlation tables, SEM coefficients, and R² summaries

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
# Load required packages
library(dplyr)
library(tidyr)
library(purrr)
library(Hmisc)
library(piecewiseSEM)
library(nlme)
library(ggplot2)
library(scales)


source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()


root        <- here::here()
data_dir    <- file.path(root, "Data", "9_SEM")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------
# 1. Import data
# ---------------------------------------------------------------------
sem_dat <- read.csv(file.path(data_dir, "SEM-data_all.csv"))  # contains Moisture_log1p, Other_biomass, etc.
str(sem_dat)

# ---------------------------------------------------------------------
# 2. Correlation structure (by species)
# ---------------------------------------------------------------------

# Variables to include in correlations
corr_vars <- c("Moisture", "SLA", "LDMC",
               "Richness_Bac_Root", "Richness_Bac_Rhizo",
               "Richness_Fun_Leaf", "Richness_Fun_Root", "Richness_Fun_Rhizo",
               "PC1_Bac_Root", "PC1_Bac_Rhizo",
               "PC1_Fun_Leaf", "PC1_Fun_Root", "PC1_Fun_Rhizo",
               "PC2_Bac_Root", "PC2_Bac_Rhizo",
               "PC2_Fun_Leaf", "PC2_Fun_Root", "PC2_Fun_Rhizo", 
               "Other_biomass", "Species_biomass")

# Helper: flatten correlation + p-value matrices
flatten_corr_matrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(Var1    = rownames(cormat)[row(cormat)[ut]],
             Var2    = rownames(cormat)[col(cormat)[ut]],
             r_value = cormat[ut], p_value = pmat[ut], row.names = NULL)
}

# Helper: compute correlations for one species
compute_species_corr <- function(df, species_name, vars) {
  sub <- df %>%
    filter(Species == species_name) %>%
    select(any_of(vars)) %>%
    select(where(is.numeric))
  
  if (nrow(sub) < 3 || ncol(sub) < 2) {
    warning("Not enough data for species: ", species_name)
    return(NULL)
  }
  
  cor_sp <- rcorr(as.matrix(sub), type = "pearson")
  out    <- flatten_corr_matrix(cor_sp$r, cor_sp$P)
  out$SpeciesGroup <- species_name
  out
}

res_festuca <- compute_species_corr(sem_dat, "Festuca", corr_vars)
res_lolium  <- compute_species_corr(sem_dat, "Lolium",  corr_vars)

# Combine and keep ordering of variables that are present
present_vars <- sort(unique(c(res_festuca$Var1, res_festuca$Var2,
                              res_lolium$Var1,  res_lolium$Var2)))
plot_vars <- corr_vars[corr_vars %in% present_vars]

alpha <- 0.05

res_all <- bind_rows(res_festuca, res_lolium) %>%
  mutate(Var1 = factor(Var1, levels = plot_vars),
         Var2 = factor(Var2, levels = plot_vars)) %>%
  filter(!is.na(Var1), !is.na(Var2)) %>%
  filter(as.integer(Var1) < as.integer(Var2))

# ---------------------------------------------------------------------
# 3. Correlation heatmap (significant pairs only)
# ---------------------------------------------------------------------

# Theme
theme_corr <- theme_classic() +
  theme(axis.text = element_text(face = "bold", color = "black", size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(face = "bold.italic", hjust = 0.5, size = 12),
        strip.background = element_blank(),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 9),
        legend.position = "inside",
        legend.position.inside = c(0.92, 0.2),
        panel.border = element_rect(fill = "black", linetype = "dotted", linewidth = 0.6))

# Build a complete upper-triangular grid for each species
grid_all <- expand_grid(Var1 = factor(plot_vars, levels = plot_vars),
                        Var2 = factor(plot_vars, levels = plot_vars)) %>%
  filter(as.integer(Var1) < as.integer(Var2))

all_half <- res_all %>%
  right_join(grid_all %>% crossing(SpeciesGroup = unique(res_all$SpeciesGroup)),
             by = c("Var1", "Var2", "SpeciesGroup")) %>%
  mutate(sig   = !is.na(p_value) & p_value < alpha,
         r_lab = ifelse(sig, sprintf("%.2f", r_value), ""))

p_corr <- ggplot(all_half, aes(Var1, Var2)) +
  # outlines for all cells
  geom_tile(width = 0.95, height = 0.95, fill = NA, color = "gray60", linewidth = 0.3) +
  # fill only significant cells
  geom_tile(data = filter(all_half, sig), aes(fill = r_value), width = 0.95, height = 0.95) +
  # labels on significant cells
  geom_text(data = filter(all_half, sig), aes(label = r_lab), size = 3, fontface = "bold") +
  scale_fill_gradient2(limits = c(-1, 1), oob = squish, name = "Pearson r",
                       low = "#0571b0", mid = "#F7F7F7", high = "#B2182B") +
  coord_equal() +
  facet_wrap(~ SpeciesGroup) +
  labs(x = NULL, y = NULL) +
  theme_corr

p_corr

ggsave(file.path(results_dir, "Significant_Correlations_bySpecies.png"), ## Figure S5
       p_corr, width = 12, height = 7, units = "in", bg = "white", dpi = 600)

# ---------------------------------------------------------------------
# 4. Export list of significant correlations (by species)
# ---------------------------------------------------------------------

sig_results_species <- sem_dat %>%
  group_by(Species) %>%
  group_modify(~ {
    sub <- .x %>%
      select(any_of(corr_vars)) %>%
      select(where(is.numeric))
    
    if (nrow(sub) < 3 || ncol(sub) < 2) return(NULL)
    
    cor_test <- rcorr(as.matrix(sub), type = "pearson")
    
    r_df <- as.data.frame(as.table(cor_test$r))
    p_df <- as.data.frame(as.table(cor_test$P))
    
    merged <- left_join(r_df, p_df, by = c("Var1", "Var2"), suffix = c("_r", "_p")) %>%
      mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
      filter(Var1 != Var2) %>%
      mutate(Vmin = pmin(Var1, Var2), Vmax = pmax(Var1, Var2)) %>%
      distinct(Vmin, Vmax, .keep_all = TRUE) %>%
      filter(Freq_p < 0.05) %>%
      transmute(Var1 = Vmin, Var2 = Vmax, r_value = round(Freq_r, 2), p_value = signif(Freq_p, 3))
    
    merged
  }) %>%
  ungroup() %>%
  mutate(SpeciesGroup = Species) %>%
  select(SpeciesGroup, Var1, Var2, r_value, p_value)

sig_results_species

write.csv(sig_results_species, row.names = FALSE,
          file.path(results_dir, "Significant_Correlations_bySpecies_raw.csv"))


# ---------------------------------------------------------------------
# 5. Helper to prepare species-specific SEM data
# ---------------------------------------------------------------------

prepare_sem_data <- function(df, species_name) {
  sub <- df %>%
    filter(Species == species_name) %>%
    mutate(Mesocosm   = factor(Mesocosm),
           Climate    = factor(Climate),
           Climate_bin = ifelse(Climate == "Future", 1, 0),
           Moisture_z  = as.numeric(scale(Moisture)))
  sub
}

# ---------------------------------------------------------------------
# 6. SEM: Festuca
# ---------------------------------------------------------------------

sub_f <- prepare_sem_data(sem_dat, "Festuca")

# Moisture_z × Climate_bin → traits / microbial variables
m_f_LDMC      <- nlme::lme(LDMC ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")
m_f_OtherBiom <- stats::lm(Other_biomass ~ Moisture_z * Climate_bin,
                           data = sub_f)
m_f_RB_Root   <- nlme::lme(Richness_Bac_Root ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")
m_f_PC1_FR    <- nlme::lme(PC1_Fun_Root ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")
m_f_PC1_FSh   <- nlme::lme(PC1_Fun_Rhizo ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")
m_f_PC2_BR    <- nlme::lme(PC2_Bac_Root ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")
m_f_PC2_BSh   <- nlme::lme(PC2_Bac_Rhizo ~ Moisture_z * Climate_bin,
                           random = ~1 | Mesocosm, data = sub_f, method = "ML")

# Biomass model
m_f_Biomass <- stats::lm(Species_biomass ~ Moisture_z * Climate_bin +
                           PC1_Bac_Rhizo + PC2_Bac_Root + Other_biomass,
                         data = sub_f)

sem_f <- piecewiseSEM::psem(m_f_LDMC, m_f_OtherBiom, m_f_RB_Root,
                            m_f_PC1_FR, m_f_PC1_FSh, m_f_PC2_BR, m_f_PC2_BSh,
                            m_f_Biomass,
                            # correlated residuals:
                            PC1_Fun_Rhizo %~~% PC1_Bac_Rhizo,
                            PC2_Bac_Rhizo %~~% PC1_Fun_Rhizo,
                            PC2_Bac_Rhizo %~~% PC2_Bac_Root,
                            PC1_Fun_Root  %~~% PC1_Bac_Rhizo,
                            PC2_Bac_Root  %~~% PC1_Fun_Rhizo,
                            data = sub_f)

summary(sem_f)
plot(sem_f)

# ---------------------------------------------------------------------
# 7. SEM: Lolium
# ---------------------------------------------------------------------

sub_l <- prepare_sem_data(sem_dat, "Lolium")

# Moisture_z × Climate_bin → traits / microbial variables
m_l_RF_Sh  <- nlme::lme(Richness_Fun_Rhizo ~ Moisture_z * Climate_bin,
                        random = ~1 | Mesocosm, data = sub_l, method = "ML")
m_l_RF_R   <- nlme::lme(Richness_Fun_Root ~ Moisture_z * Climate_bin,
                        random = ~1 | Mesocosm, data = sub_l, method = "ML")
m_l_Other  <- stats::lm(Other_biomass ~ Moisture_z * Climate_bin,
                        data = sub_l)
m_l_PC1_FR <- nlme::lme(PC1_Fun_Root ~ Moisture_z * Climate_bin,
                        random = ~1 | Mesocosm, data = sub_l, method = "ML")

# Biomass model
m_l_Biomass <- stats::lm(Species_biomass ~ Moisture_z * Climate_bin +
                           Other_biomass + SLA, data = sub_l)

sem_l <- piecewiseSEM::psem(m_l_RF_Sh, m_l_RF_R, m_l_Other, m_l_PC1_FR, m_l_Biomass,
                            # correlated residuals:
                            Other_biomass     %~~% SLA,
                            Richness_Fun_Root %~~% Richness_Fun_Rhizo,
                            data = sub_l)

summary(sem_l)
plot(sem_l)

# ---------------------------------------------------------------------
# 8. Extract and save standardized coefficients and R²
# ---------------------------------------------------------------------

# Standardized path coefficients
coefs_festuca <- coefs(sem_f, standardize = "scale")
coefs_lolium  <- coefs(sem_l, standardize = "scale")

coefs_festuca$Species <- "Festuca"
coefs_lolium$Species  <- "Lolium"

coefs_all <- rbind(coefs_festuca, coefs_lolium)

write.csv(coefs_all, row.names = FALSE,
          file.path(results_dir, "SEM_Coefs_Comparison.csv")) # Table S14


# Marginal / conditional R² for component models
r2_festuca <- rsquared(sem_f)
r2_lolium  <- rsquared(sem_l)

r2_festuca$Species <- "Festuca"
r2_lolium$Species  <- "Lolium"

r2_all <- rbind(r2_festuca, r2_lolium)

write.csv(r2_all, row.names = FALSE,
          file.path(results_dir, "SEM_R2_Comparison.csv"))     # Table S15

