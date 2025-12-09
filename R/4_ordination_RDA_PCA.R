##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/4_ordination_RDA_PCA.R
##
## Tasks:
##  1) Hellinger-transform ASV tables for 6 compartments (BL, BR, BS, FL, FR, FS)
##  2) RDA with Moisture × Climate × Species + permutation tests
##       - variance explained (%) per term (VarExp_RDA.csv)
##       - term-level ANOVA tables for Moisture (RDA_terms_Moisture.csv)
##  3) RDA with GWC × Climate × Species + permutation tests
##       - term-level ANOVA tables for GWC (RDA_terms_GWC.csv)
##  4) Unconstrained ordination (PCA via rda(~1)) for each compartment
##       - PCA plots colored by Plant species (PCA_all_bySpecies.png)
##       - PCA plots colored by Climate (PCA_all_byClimate.png)

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
ps_dir      <- file.path(root, "Data", "2_Phyloseq objects", "Output", "ps_new")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Ordination theme
theme_ord <- theme_classic() +
  theme(legend.spacing.x = unit(0.01, "lines"),
        axis.line        = element_line(color = "black"),
        title            = element_text(face = "bold", size = 12),
        axis.title       = element_text(face = "bold", size = 12),
        legend.title     = element_text(face = "bold", size = 12),
        legend.text      = element_text(face = "bold", size = 12),
        plot.title       = element_text(face = "bold", size = 12))

# Convenience palettes from utils_theme.R
species_cols <- species_cols
climate_cols <- climate_cols

# -------------------------------------------------------------------
# 1) Load phyloseq objects and prepare ASV/env tables
# -------------------------------------------------------------------
ps_BL <- readRDS(file.path(ps_dir, "ps-BL.rds"))  # Bacteria - Leaf
ps_BR <- readRDS(file.path(ps_dir, "ps-BR.rds"))  # Bacteria - Root
ps_BS <- readRDS(file.path(ps_dir, "ps-BS.rds"))  # Bacteria - Rhizosphere
ps_FL <- readRDS(file.path(ps_dir, "ps-FL.rds"))  # Fungi    - Leaf
ps_FR <- readRDS(file.path(ps_dir, "ps-FR.rds"))  # Fungi    - Root
ps_FS <- readRDS(file.path(ps_dir, "ps-FS.rds"))  # Fungi    - Rhizosphere

# helper: extract ASV matrix (samples x taxa) and env data.frame
extract_asv_env <- function(ps) {
  asv <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) asv <- t(asv)
  env <- as.data.frame(phyloseq::sample_data(ps))
  
  env$Mesocosm <- factor(env$Mesocosm)
  env$Moisture <- as.numeric(env$Moisture)
  if ("GWC" %in% names(env)) env$GWC <- as.numeric(env$GWC)
  
  list(asv = asv, env = env)
}

dat_BL <- extract_asv_env(ps_BL)
dat_BR <- extract_asv_env(ps_BR)
dat_BS <- extract_asv_env(ps_BS)
dat_FL <- extract_asv_env(ps_FL)
dat_FR <- extract_asv_env(ps_FR)
dat_FS <- extract_asv_env(ps_FS)

# Hellinger transformation
asv_BLh <- vegan::decostand(dat_BL$asv, method = "hellinger")
asv_BRh <- vegan::decostand(dat_BR$asv, method = "hellinger")
asv_BSh <- vegan::decostand(dat_BS$asv, method = "hellinger")
asv_FLh <- vegan::decostand(dat_FL$asv, method = "hellinger")
asv_FRh <- vegan::decostand(dat_FR$asv, method = "hellinger")
asv_FSh <- vegan::decostand(dat_FS$asv, method = "hellinger")

env1 <- dat_BL$env
env2 <- dat_BR$env
env3 <- dat_BS$env
env4 <- dat_FL$env
env5 <- dat_FR$env
env6 <- dat_FS$env

# -------------------------------------------------------------------
# 2) RDA: Moisture × Climate × Species
# -------------------------------------------------------------------
set.seed(2025)

# Permutation controls (Mesocosm as strata where appropriate)
CTRL1 <- permute::how(nperm = 999)

CTRL2 <- permute::how(nperm = 999,
                      plots  = permute::Plots(strata = env2$Mesocosm, type = "free"),
                      within = permute::Within(type = "none"))

CTRL3 <- permute::how(nperm = 999,
                      plots  = permute::Plots(strata = env3$Mesocosm, type = "free"),
                      within = permute::Within(type = "none"))

CTRL4 <- permute::how(nperm = 999)

CTRL5 <- permute::how(nperm = 999,
                      plots  = permute::Plots(strata = env5$Mesocosm, type = "free"),
                      within = permute::Within(type = "none"))

CTRL6 <- permute::how(nperm = 999,
                      plots  = permute::Plots(strata = env6$Mesocosm, type = "free"),
                      within = permute::Within(type = "none"))

# RDA models (Moisture only)
b1 <- vegan::rda(asv_BLh ~ Moisture * Climate * Species, data = env1)
b2 <- vegan::rda(asv_BRh ~ Moisture * Climate * Species, data = env2)
b3 <- vegan::rda(asv_BSh ~ Moisture * Climate * Species, data = env3)
b4 <- vegan::rda(asv_FLh ~ Moisture * Climate * Species, data = env4)
b5 <- vegan::rda(asv_FRh ~ Moisture * Climate * Species, data = env5)
b6 <- vegan::rda(asv_FSh ~ Moisture * Climate * Species, data = env6)

# ANOVA by terms
ab1 <- vegan::anova.cca(b1, by = "terms", permutations = CTRL1)
ab2 <- vegan::anova.cca(b2, by = "terms", permutations = CTRL2)
ab3 <- vegan::anova.cca(b3, by = "terms", permutations = CTRL3)
ab4 <- vegan::anova.cca(b4, by = "terms", permutations = CTRL4)
ab5 <- vegan::anova.cca(b5, by = "terms", permutations = CTRL5)
ab6 <- vegan::anova.cca(b6, by = "terms", permutations = CTRL6)

rda_models_moist <- list(BL = b1, BR = b2, BS = b3, FL = b4, FR = b5, FS = b6)
rda_terms_moist  <- list(BL = ab1, BR = ab2, BS = ab3, FL = ab4, FR = ab5, FS = ab6)

# -------------------------------------------------------------------
# 3) Variance explained (%) per term (Moisture models)
# -------------------------------------------------------------------
compute_rda_var_exp <- function(mod, aov_tab, comp_label) {
  df <- as.data.frame(aov_tab)
  df$Term <- rownames(df)
  rownames(df) <- NULL
  
  df <- df[df$Term != "Residual", , drop = FALSE]
  if (!"Variance" %in% names(df)) {
    stop("ANOVA table does not contain a 'Variance' column.")
  }
  
  total_var <- sum(df$Variance, na.rm = TRUE)
  R2adj     <- vegan::RsquareAdj(mod)$adj.r.squared
  
  df$Compartment <- comp_label
  df$Variance_explained_percent <- 100 * (df$Variance / total_var) * R2adj
  df$Variance_explained_percent <- round(df$Variance_explained_percent, 4)
  df$R2adj <- round(R2adj, 4)
  
  df[, c("Compartment", "Term", "Variance_explained_percent", "R2adj")]
}

var_exp_rda_list <- mapply(FUN  = compute_rda_var_exp,
                           mod  = rda_models_moist,
                           aov_tab = rda_terms_moist,
                           comp = names(rda_models_moist),
                           SIMPLIFY = FALSE)

var_exp_rda_tbl <- do.call(rbind, var_exp_rda_list)

write.csv(var_exp_rda_tbl, file.path(results_dir, "VarExp_RDA.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4) Term-level ANOVA tables: Moisture & GWC
# -------------------------------------------------------------------

# Helper to tidy ANOVA tables
tidy_anova <- function(aov_tab, dataset, driver) {
  df <- as.data.frame(aov_tab)
  df$term <- rownames(df)
  rownames(df) <- NULL
  
  cols_to_keep <- intersect(c("Df", "Variance", "F", "Pr(>F)"), names(df))
  df_out <- df[, c("term", cols_to_keep), drop = FALSE]
  df_out$dataset <- dataset
  df_out$driver  <- driver
  
  df_out[, c("dataset", "driver", "term", cols_to_keep)]
}

# Moisture term tables
moisture_results <- rbind(
  tidy_anova(ab1, "BL (Bacteria-Leaf)",  "Moisture"),
  tidy_anova(ab2, "BR (Bacteria-Root)",  "Moisture"),
  tidy_anova(ab3, "BS (Bacteria-Soil)",  "Moisture"),
  tidy_anova(ab4, "FL (Fungi-Leaf)",     "Moisture"),
  tidy_anova(ab5, "FR (Fungi-Root)",     "Moisture"),
  tidy_anova(ab6, "FS (Fungi-Soil)",     "Moisture")
)

write.csv(moisture_results, file.path(results_dir, "RDA_terms_Moisture.csv"), ## Table 3
          row.names = FALSE)

# -------------------------------------------------------------------
# 5) RDA: GWC × Climate × Species
# -------------------------------------------------------------------
set.seed(2025)

c1 <- vegan::rda(asv_BLh ~ GWC * Climate * Species, data = env1)
c2 <- vegan::rda(asv_BRh ~ GWC * Climate * Species, data = env2)
c3 <- vegan::rda(asv_BSh ~ GWC * Climate * Species, data = env3)
c4 <- vegan::rda(asv_FLh ~ GWC * Climate * Species, data = env4)
c5 <- vegan::rda(asv_FRh ~ GWC * Climate * Species, data = env5)
c6 <- vegan::rda(asv_FSh ~ GWC * Climate * Species, data = env6)

ac1 <- vegan::anova.cca(c1, by = "terms", permutations = CTRL1)
ac2 <- vegan::anova.cca(c2, by = "terms", permutations = CTRL2)
ac3 <- vegan::anova.cca(c3, by = "terms", permutations = CTRL3)
ac4 <- vegan::anova.cca(c4, by = "terms", permutations = CTRL4)
ac5 <- vegan::anova.cca(c5, by = "terms", permutations = CTRL5)
ac6 <- vegan::anova.cca(c6, by = "terms", permutations = CTRL6)

gwc_results <- rbind(
  tidy_anova(ac1, "BL (Bacteria-Leaf)",  "GWC"),
  tidy_anova(ac2, "BR (Bacteria-Root)",  "GWC"),
  tidy_anova(ac3, "BS (Bacteria-Soil)",  "GWC"),
  tidy_anova(ac4, "FL (Fungi-Leaf)",     "GWC"),
  tidy_anova(ac5, "FR (Fungi-Root)",     "GWC"),
  tidy_anova(ac6, "FS (Fungi-Soil)",     "GWC")
)

write.csv(gwc_results, file.path(results_dir, "RDA_terms_GWC.csv"),  ## Table S7
          row.names = FALSE)

# -------------------------------------------------------------------
# 6) PCA (unconstrained) via RDA ~ 1
# -------------------------------------------------------------------
a1 <- vegan::rda(asv_BLh)  # Bacteria - Leaf
a2 <- vegan::rda(asv_BRh)  # Bacteria - Root
a3 <- vegan::rda(asv_BSh)  # Bacteria - Rhizosphere
a4 <- vegan::rda(asv_FLh)  # Fungi    - Leaf
a5 <- vegan::rda(asv_FRh)  # Fungi    - Root
a6 <- vegan::rda(asv_FSh)  # Fungi    - Rhizosphere

# -------------------------------------------------------------------
# 7) PCA plot helpers (convex hulls, species-colored, climate-colored)
# -------------------------------------------------------------------

# Convex hull per group
hull_df <- function(dat, grp_col = "Species") {
  dat[[grp_col]] <- factor(dat[[grp_col]])
  sp <- split(seq_len(nrow(dat)), dat[[grp_col]])
  inds <- lapply(sp, function(ii) {
    pts <- dat[ii, c("X", "Y")]
    pts <- pts[stats::complete.cases(pts), , drop = FALSE]
    if (nrow(unique(pts)) >= 3) {
      ii[grDevices::chull(pts$X, pts$Y)]
    } else ii
  })
  do.call(rbind,
          lapply(names(inds), function(g) {
            dat[inds[[g]], c("X", "Y", grp_col), drop = FALSE]
            }))
}

# Generic: build site scores + env into long df
build_pca_df <- function(ord, env, scaling = 2) {
  scr   <- vegan::scores(ord, display = "sites", scaling = scaling)
  df    <- data.frame(scr, SampleID = rownames(scr))
  env_df <- as.data.frame(env)
  env_df$SampleID <- rownames(env_df)
  
  dat <- merge(df, env_df[, c("SampleID", "Climate", "Species", "Moisture")],
               by = "SampleID", all.x = TRUE, sort = FALSE)
  
  dat$Species <- factor(dat$Species)
  dat$Climate <- factor(dat$Climate)
  
  ax <- colnames(scr)[1:2]
  dat$X <- dat[[ax[1]]]
  dat$Y <- dat[[ax[2]]]
  
  imp <- summary(ord)$cont$importance[2, 1:2] * 100
  xlab_txt <- sprintf("%s (%.1f%%)", ax[1], imp[1])
  ylab_txt <- sprintf("%s (%.1f%%)", ax[2], imp[2])
  
  list(dat = dat, xlab = xlab_txt, ylab = ylab_txt)
}

# PCA: color by Species, size = Moisture
pca_plot_species <- function(ord, env, title, scaling = 2,
                             species_colors = species_cols) {
  dd <- build_pca_df(ord, env, scaling = scaling)
  dat  <- dd$dat
  xlab_txt <- dd$xlab
  ylab_txt <- dd$ylab
  
  hulls <- hull_df(dat, "Species")
  
  ggplot(dat, aes(X, Y)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_polygon(data = hulls, aes(group = Species, color = Species),
                 fill = NA, linewidth = 0.9, show.legend = FALSE) +
    geom_point(aes(color = Species, size = Moisture), alpha = 0.7) +
    scale_color_manual(values = species_colors, name = "Plant species") +
    scale_size_continuous(name = "Moisture (ml)") +
    labs(title = title, x = xlab_txt, y = ylab_txt) +
    theme_ord +
    guides(color = guide_legend(override.aes = list(size = 3),
                                order        = 1,
                                title.theme  = element_text(face = "bold"),
                                label.theme  = element_text(face = "bold.italic", size = 12)),
           size = guide_legend(title.theme = element_text(face = "bold"), order = 2))
}

# PCA: color by Climate, size = Moisture
pca_plot_climate <- function(ord, env, title, scaling = 2,
                             climate_colors = climate_cols) {
  dd <- build_pca_df(ord, env, scaling = scaling)
  dat  <- dd$dat
  xlab_txt <- dd$xlab
  ylab_txt <- dd$ylab
  
  hulls <- hull_df(dat, "Climate")
  
  ggplot(dat, aes(X, Y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_polygon(data = hulls, aes(group = Climate, color = Climate),
                 fill = NA, linewidth = 0.9, show.legend = FALSE) +
    geom_point(aes(color = Climate, size = Moisture), alpha = 0.7) +
    scale_color_manual(values = climate_colors, name = "Climatic conditions") +
    scale_size_continuous(name = "Moisture (ml)") +
    labs(title = title, x = xlab_txt, y = ylab_txt) +
    theme_ord +
    guides(color = guide_legend(override.aes = list(size = 3),
                                title.theme  = element_text(face = "bold"),
                                label.theme  = element_text(face = "bold", size = 12)),
           size = guide_legend(title.theme = element_text(face = "bold")))
}

# -------------------------------------------------------------------
# 8) Build PCA panels and save figures
# -------------------------------------------------------------------

# Species-colored PCA
p1_sp <- pca_plot_species(a1, env1, "Bacteria · Leaf endosphere (BL)")
p2_sp <- pca_plot_species(a2, env2, "Bacteria · Root endosphere (BR)")
p3_sp <- pca_plot_species(a3, env3, "Bacteria · Rhizosphere (BS)")
p4_sp <- pca_plot_species(a4, env4, "Fungi · Leaf endosphere (FL)")
p5_sp <- pca_plot_species(a5, env5, "Fungi · Root endosphere (FR)")
p6_sp <- pca_plot_species(a6, env6, "Fungi · Rhizosphere (FS)")

p_all_sp <- ggpubr::ggarrange(p1_sp, p2_sp, p3_sp, p4_sp, p5_sp, p6_sp,
                              nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom",
                              labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                              font.label  = list(size = 12))

ggplot2::ggsave(file.path(results_dir, "PCA_all_bySpecies.png"), ## Figure 3
                p_all_sp, width = 14, height = 8, units = "in", bg = "white", dpi = 600)

# Climate-colored PCA
p1_cl <- pca_plot_climate(a1, env1, "Bacteria · Leaf endosphere (BL)")
p2_cl <- pca_plot_climate(a2, env2, "Bacteria · Root endosphere (BR)")
p3_cl <- pca_plot_climate(a3, env3, "Bacteria · Rhizosphere (BS)")
p4_cl <- pca_plot_climate(a4, env4, "Fungi · Leaf endosphere (FL)")
p5_cl <- pca_plot_climate(a5, env5, "Fungi · Root endosphere (FR)")
p6_cl <- pca_plot_climate(a6, env6, "Fungi · Rhizosphere (FS)")

p_all_cl <- ggpubr::ggarrange(p1_cl, p2_cl, p3_cl, p4_cl, p5_cl, p6_cl,
                              nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom",
                              labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                              font.label  = list(size = 12))

ggplot2::ggsave(file.path(results_dir, "PCA_all_byClimate.png"), ## Figure S10
                p_all_cl, width = 14, height = 8, units = "in", bg = "white", dpi = 600)

