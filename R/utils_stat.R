##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## R/utils_stats.R

# Standardize continuous predictor and ensure factors exist
prep_data_z <- function(dat, cont_var, random_factor = NULL,
                        factor_vars = c("Climate", "Species", "Mesocosm")) {
  for (f in factor_vars) {
    if (f %in% names(dat)) {
      dat[[f]] <- factor(dat[[f]])
    }
  }
  zname <- paste0(cont_var, "_z")
  dat[[zname]] <- scale(as.numeric(dat[[cont_var]]), center = TRUE, scale = TRUE)
  dat
}

# Fit a generic LMM with response ~ cont * Climate * Species + (1 | random_factor)
fit_lmm_3way <- function(dat, response, cont_var, random_factor = "Mesocosm") {
  dat <- prep_data_z(dat, cont_var = cont_var, random_factor = random_factor)
  zname <- paste0(cont_var, "_z")
  
  form <- as.formula(sprintf("%s ~ %s * Climate * Species + (1 | %s)",
                             response, zname, random_factor))
  lme4::lmer(form, data = dat)
}

# Type-III ANOVA for a fitted model
anova_type3 <- function(mod) {
  car::Anova(mod, type = "III")
}

# Helper for standardizing Anova tables
std_names_anova <- function(df) {
  if ("F value" %in% names(df)) names(df)[names(df) == "F value"] <- "Statistic"
  if ("Chisq"   %in% names(df)) names(df)[names(df) == "Chisq"]   <- "Statistic"
  if ("Pr(>F)" %in% names(df))  names(df)[names(df) == "Pr(>F)"]  <- "p_value"
  if ("Pr(>Chisq)" %in% names(df)) names(df)[names(df) == "Pr(>Chisq)"] <- "p_value"
  df
}

# Variance explained (%) for LMM Chi-square tables (fixed effects)
variance_explained_lmm <- function(anova_tab, model) {
  df <- as.data.frame(anova_tab)
  df$Term <- rownames(df)
  df <- df[df$Term != "(Intercept)", , drop = FALSE]
  
  chisq_total <- sum(df[["Chisq"]], na.rm = TRUE)
  R2m <- MuMIn::r.squaredGLMM(model)[, "R2m"]
  
  df$Variance_explained_percent <- 100 * (df[["Chisq"]] / chisq_total) * R2m
  df$Variance_explained_percent <- round(df$Variance_explained_percent, 4)
  
  df$Response <- deparse(formula(model)[[2]])
  df[, c("Response", "Term", "Variance_explained_percent")]
}
