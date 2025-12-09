##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/1b_plant_performance_visualization.R
##
## Visualization of plant performance:
## - Moisture gradients with drought zones
## - Species / Climate / Species×Climate violins
##
## Inputs:
##   Data/1_Plant performance/Species biomass.csv
##   Data/1_Plant performance/SLA+LDMC.csv
##
## Outputs (saved in ./figures/):
##   Plant_performance_Species.png
##   Plant_performance_MxC.png
##   Plant_performance_Climate.png

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root       <- here::here()
data_dir   <- file.path(root, "Data", "1_Plant performance")
results_dir <- file.path(root, "results")
fig_dir     <- file.path(root, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)

theme_facet <- project_theme_facet()

cs <- species_cols   # from utils_theme.R
cc <- climate_cols   # from utils_theme.R

# -------------------------------------------------------------------
# 1) Load and prepare data
# -------------------------------------------------------------------
data_biomass <- read.csv(file.path(data_dir, "Species biomass.csv"), header = TRUE)
data_traits  <- read.csv(file.path(data_dir, "SLA+LDMC.csv"),        header = TRUE)

data_biomass <- data_biomass |>
  dplyr::mutate(Species  = factor(Species),
                Climate  = factor(Climate, levels = c("Current", "Future")),
                Mesocosm = factor(Mesocosm))

data_traits <- data_traits |>
  dplyr::mutate(Species  = factor(Species),
                Climate  = factor(Climate, levels = c("Current", "Future")),
                Mesocosm = factor(Mesocosm))

species_levels <- levels(data_biomass$Species)
climate_levels <- levels(data_biomass$Climate)

# -------------------------------------------------------------------
# 2) Long data for Moisture-gradient plots (4 responses)
# -------------------------------------------------------------------
bio_biomass <- data_biomass |>
  dplyr::transmute(Moisture, Species, Climate,
                   Response = "Species biomass",
                   Value    = Species_biomass)

bio_prop <- data_biomass |>
  dplyr::transmute(Moisture, Species, Climate,
                   Response = "Biomass proportion",
                   Value    = Proportion)

trait_sla <- data_traits |>
  dplyr::transmute(Moisture, Species, Climate,
                   Response = "SLA",
                   Value    = SLA)

trait_ldmc <- data_traits |>
  dplyr::transmute(Moisture, Species, Climate,
                   Response = "LDMC",
                   Value    = LDMC)

plot_dat <- dplyr::bind_rows(bio_biomass, bio_prop, trait_sla, trait_ldmc) |>
  dplyr::mutate(Response = factor(
      Response, levels = c("Species biomass", "Biomass proportion", "SLA", "LDMC")))

# Moisture range + drought thresholds
x_min <- min(plot_dat$Moisture, na.rm = TRUE)
x_max <- max(plot_dat$Moisture, na.rm = TRUE)
thr1  <- 400
thr2  <- 1000

# Top strip with drought labels (different geometry than add_drought_zones)
make_drought_top_strip <- function(x_min, x_max, thr1, thr2,
                                   labels = c("Severe drought",
                                              "Moderate drought",
                                              "Well-watered")) {
  mid_severe   <- (x_min + thr1) / 2
  mid_moderate <- (thr1 + thr2) / 2
  mid_wet      <- (thr2 + x_max) / 2
  
  lab_df <- data.frame(x     = c(mid_severe, mid_moderate, mid_wet),
                       y     = 0.5, label = labels)
  
  ggplot(lab_df, aes(x = x, y = y, label = label)) +
    # same colors as add_drought_zones, but spanning 0–1
    annotate("rect", xmin = x_min, xmax = thr1, ymin = 0, ymax = 1,
             fill = "#f9d7c2", alpha = 0.4) +
    annotate("rect", xmin = thr1, xmax = thr2, ymin = 0, ymax = 1,
             fill = "#fffec0", alpha = 0.4) +
    annotate("rect", xmin = thr2, xmax = x_max, ymin = 0, ymax = 1,
             fill = "#d9eaf7", alpha = 0.4) +
    geom_text(fontface = "bold", size = 3.5) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 1, r = 15, b = -10, l = 15))
}

# -------------------------------------------------------------------
# 3) Moisture × Species (4 responses, pooled across climate)
# -------------------------------------------------------------------

## 3.1 Slopes per Response × Species (no broom, do it by hand)
slopes_species <- plot_dat |>
  dplyr::group_by(Response, Species) |>
  dplyr::group_modify(~ {
    m <- lm(Value ~ Moisture, data = .x)
    sm <- summary(m)$coefficients
    est <- sm["Moisture", "Estimate"]
    p   <- sm["Moisture", "Pr(>|t|)"]
    dplyr::tibble(estimate = est, p.value = p)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                stars     = p_to_stars(p.value),
                label     = paste0("Slope = ",
                                   formatC(estimate, format = "f", digits = 3),
                                   " ", stars))

plot_dat_species <- plot_dat |>
  dplyr::left_join(slopes_species |>
                     dplyr::select(Response, Species, line_type, label),
                   by = c("Response", "Species"))

panel_range <- plot_dat_species |>
  dplyr::group_by(Response) |>
  dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

annot_slopes_species <- slopes_species |>
  dplyr::left_join(panel_range, by = "Response") |>
  dplyr::group_by(Response) |>
  dplyr::arrange(Species) |>
  dplyr::mutate(x = x_max * 0.64,
                y = y_max * (0.95 - 0.10 * (dplyr::row_number() - 1))) |>
  dplyr::ungroup()

## 3.2 Main Moisture × Species plot (4 rows)
p_raw_sig <- ggplot(plot_dat_species, aes(x = Moisture, y = Value, color = Species)) +
  add_drought_zones(x_min, x_max, thr1, thr2) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE,
              aes(group = Species, linetype = line_type, fill = Species),
              alpha = 0.15, color = NA, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE,
              aes(group = Species, linetype = line_type, color = Species),
              linewidth = 1.0) +
  facet_grid(Response ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_color_manual(values = cs) +
  scale_fill_manual(values = cs) +
  scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed")) +
  guides(linetype = guide_legend(title = "Significance", order = 2,
                                 override.aes = list(color = "black",
                                                     linetype = c("dashed", "solid"),
                                                     linewidth = 0.5)),
    color = guide_legend(title = "Plant species", order = 1)) +
  labs(x = "Moisture (watering amount, ml)", y = "Observed value") +
  theme_facet +
  geom_text(data = annot_slopes_species,
            aes(x = x, y = y, label = label, color = Species),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1,
            size = 3.5, fontface = "bold", show.legend = FALSE)

p_top <- make_drought_top_strip(x_min, x_max, thr1, thr2)
final_plot_MP <- p_top / p_raw_sig + patchwork::plot_layout(heights = c(0.06, 1))


# -------------------------------------------------------------------
# 4) Moisture-only (no Species/Climate) – slopes per Response
# -------------------------------------------------------------------
slopes_M <- plot_dat |>
  dplyr::group_by(Response) |>
  dplyr::group_modify(~ {
    m <- lm(Value ~ Moisture, data = .x)
    sm <- summary(m)$coefficients
    est <- sm["Moisture", "Estimate"]
    p   <- sm["Moisture", "Pr(>|t|)"]
    dplyr::tibble(estimate = est, p.value = p)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(stars = p_to_stars(p.value),
                label = paste0("Slope = ",
                               formatC(estimate, format = "f", digits = 3),
                               " ", stars))

panel_range_M <- plot_dat |>
  dplyr::group_by(Response) |>
  dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

annot_M <- slopes_M |>
  dplyr::left_join(panel_range_M, by = "Response") |>
  dplyr::mutate(x = x_max * 0.64, y = y_max * 0.90)

p_M_only <- ggplot(plot_dat, aes(x = Moisture, y = Value)) +
  add_drought_zones(x_min, x_max, thr1, thr2) +
  geom_point(alpha = 0.5, size = 1.5, color = "#6A6A6A") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "grey50",
              alpha = 0.2, linewidth = 1.0, linetype = "dashed") +
  facet_grid(Response ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  labs(x = "Moisture (watering amount, ml)", y = "Observed value") +
  theme_facet +
  geom_text(data = annot_M, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5, fontface = "bold")

p_top_M <- make_drought_top_strip(x_min, x_max, thr1, thr2)
final_plot_M <- p_top_M / p_M_only + patchwork::plot_layout(heights = c(0.06, 1))


# -------------------------------------------------------------------
# 5) Moisture × Species × Climate (M × C × P)
# -------------------------------------------------------------------
slopes_RSPC <- plot_dat |>
  dplyr::group_by(Response, Species, Climate) |>
  dplyr::group_modify(~ {
    m <- lm(Value ~ Moisture, data = .x)
    sm <- summary(m)$coefficients
    est <- sm["Moisture", "Estimate"]
    p   <- sm["Moisture", "Pr(>|t|)"]
    dplyr::tibble(estimate = est, p.value = p)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                stars     = p_to_stars(p.value),
                label     = paste0("Slope = ",
                                   formatC(estimate, format = "f", digits = 3),
                                   " ", stars))

plot_dat_RSPC <- plot_dat |>
  dplyr::left_join(slopes_RSPC |>
                     dplyr::select(Response, Species, Climate, line_type, label),
                   by = c("Response", "Species", "Climate"))

panel_range_R <- plot_dat_RSPC |>
  dplyr::group_by(Response) |>
  dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

annot_RSPC <- slopes_RSPC |>
  dplyr::left_join(panel_range_R, by = "Response") |>
  dplyr::group_by(Response, Climate) |>
  dplyr::arrange(Species) |>
  dplyr::mutate(x = x_max * 0.54,
                y = y_max * (0.95 - 0.10 * (dplyr::row_number() - 1))) |>
  dplyr::ungroup()

p_raw_sig_clim <- ggplot(plot_dat_RSPC,
                         aes(x = Moisture, y = Value, color = Species)) +
  add_drought_zones(x_min, x_max, thr1, thr2) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE,
              aes(group = Species, linetype = line_type, fill = Species),
              alpha = 0.15, color = NA, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE,
              aes(group = Species, linetype = line_type, color = Species),
              linewidth = 1.0) +
  facet_grid(Response ~ Climate, scales = "free_y") +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_color_manual(values = cs) +
  scale_fill_manual(values = cs) +
  scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed")) +
  guides(linetype = guide_legend(title = "Significance", order = 2,
                                 override.aes = list(color = "black",
                                                     linetype = c("dashed", "solid"),
                                                     linewidth = 0.5)),
    color   = guide_legend(title = "Plant species", order = 1)) +
  labs(x = "Moisture (watering amount, ml)", y = "Observed value") +
  theme_facet +
  theme(legend.position = "bottom") +
  geom_text(data = annot_RSPC, aes(x = x, y = y, label = label, color = Species),
            inherit.aes = FALSE, hjust = 0, vjust = 0,
            size = 3.5, fontface = "bold", show.legend = FALSE)

p_top_MC <- make_drought_top_strip(x_min, x_max, thr1, thr2)
final_plot_MC <- p_top_MC / p_raw_sig_clim +
  patchwork::plot_layout(heights = c(0.06, 1))


# -------------------------------------------------------------------
# 6) Violin plots
#    - Species only
#    - Climate only
#    - Species × Climate
# -------------------------------------------------------------------
## 6.1 Species-only violin plots (4 responses)
bio_biomass_sp <- data_biomass |>
  dplyr::select(Mesocosm, Species, Species_biomass) |>
  dplyr::rename(Value = Species_biomass) |>
  dplyr::mutate(Response = "Species biomass")

bio_prop_sp <- data_biomass |>
  dplyr::select(Mesocosm, Species, Proportion) |>
  dplyr::rename(Value = Proportion) |>
  dplyr::mutate(Response = "Biomass proportion")

trait_sla_sp <- data_traits |>
  dplyr::select(Mesocosm, Species, SLA) |>
  dplyr::rename(Value = SLA) |>
  dplyr::mutate(Response = "SLA")

trait_ldmc_sp <- data_traits |>
  dplyr::select(Mesocosm, Species, LDMC) |>
  dplyr::rename(Value = LDMC) |>
  dplyr::mutate(Response = "LDMC")

df_violin_sp <- dplyr::bind_rows(bio_biomass_sp, bio_prop_sp,
                                 trait_sla_sp, trait_ldmc_sp) |>
  dplyr::filter(!is.na(Value)) |>
  dplyr::mutate(Species  = factor(Species),
                Response = factor(Response,
                                  levels = c("Species biomass", "Biomass proportion",
                                             "SLA", "LDMC")))

iqr_sp <- df_violin_sp |>
  dplyr::group_by(Response, Species) |>
  dplyr::summarise(median = median(Value, na.rm = TRUE),
                   q1     = quantile(Value, 0.25, na.rm = TRUE),
                   q3     = quantile(Value, 0.75, na.rm = TRUE),
                   .groups = "drop")

letters_sp <- df_violin_sp |>
  dplyr::group_by(Response) |>
  dplyr::group_modify(~ {
    resp <- .y$Response
    if (resp %in% c("SLA", "LDMC")) {
      m <- lme4::lmer(Value ~ Species + (1 | Mesocosm), data = .x)
    } else {
      m <- stats::lm(Value ~ Species, data = .x)
    }
    emm <- emmeans::emmeans(m, ~ Species)
    cl  <- multcomp::cld(emm, adjust = "none", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Species, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    df_violin_sp |>
      dplyr::group_by(Response) |>
      dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop"),
    by = "Response") |>
  dplyr::mutate(y_pos = y_max * 1.05)

p_violin_sp <- ggplot(df_violin_sp, aes(x = Species, y = Value, fill = Species)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_point(data = iqr_sp, aes(x = Species, y = median), inherit.aes = FALSE,
             size = 2, shape = 21, fill = "black", color = "black", stroke = 0.5) +
  geom_segment(data = iqr_sp, aes(x = Species, xend = Species, y = q1, yend = q3),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_text(data = letters_sp, aes(x = Species, y = y_pos, label = group_letters),
            inherit.aes = FALSE, fontface = "bold", size = 3) +
  facet_grid(Response ~ ., scales = "free_y") +
  scale_fill_manual(values = cs) +
  labs(x = "Plant species", y = "Observed value") +
  theme_facet +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold.italic", size = 10, color = "black"))

## 6.2 Climate-only violins
bio_biomass_cl <- data_biomass |>
  dplyr::select(Mesocosm, Climate, Species_biomass) |>
  dplyr::rename(Value = Species_biomass) |>
  dplyr::mutate(Response = "Species biomass")

bio_prop_cl <- data_biomass |>
  dplyr::select(Mesocosm, Climate, Proportion) |>
  dplyr::rename(Value = Proportion) |>
  dplyr::mutate(Response = "Biomass proportion")

trait_sla_cl <- data_traits |>
  dplyr::select(Mesocosm, Climate, SLA) |>
  dplyr::rename(Value = SLA) |>
  dplyr::mutate(Response = "SLA")

trait_ldmc_cl <- data_traits |>
  dplyr::select(Mesocosm, Climate, LDMC) |>
  dplyr::rename(Value = LDMC) |>
  dplyr::mutate(Response = "LDMC")

df_clim <- dplyr::bind_rows(bio_biomass_cl, bio_prop_cl, trait_sla_cl, trait_ldmc_cl) |>
  dplyr::filter(!is.na(Value)) |>
  dplyr::mutate(Climate  = factor(Climate, levels = c("Current", "Future")),
                Response = factor(Response,
                                  levels = c("Species biomass", "Biomass proportion",
                                             "SLA", "LDMC")))

iqr_cl <- df_clim |>
  dplyr::group_by(Response, Climate) |>
  dplyr::summarise(median = median(Value, na.rm = TRUE),
                   q1     = quantile(Value, 0.25, na.rm = TRUE),
                   q3     = quantile(Value, 0.75, na.rm = TRUE),
                   .groups = "drop")

letters_cl <- df_clim |>
  dplyr::group_by(Response) |>
  dplyr::group_modify(~ {
    resp <- .y$Response
    if (resp %in% c("SLA", "LDMC")) {
      m <- lme4::lmer(Value ~ Climate + (1 | Mesocosm), data = .x)
    } else {
      m <- stats::lm(Value ~ Climate, data = .x)
    }
    emm <- emmeans::emmeans(m, ~ Climate)
    cl  <- multcomp::cld(emm, adjust = "none", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Climate, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    df_clim |>
      dplyr::group_by(Response) |>
      dplyr::summarise(
        y_max = max(Value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Response"
  ) |>
  dplyr::mutate(y_pos = y_max * 1.05)

p_violin_cl <- ggplot(df_clim, aes(x = Climate, y = Value, fill = Climate)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_segment(data = iqr_cl, aes(x = Climate, xend = Climate, y = q1, yend = q3),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_point(data = iqr_cl, aes(x = Climate, y = median), inherit.aes = FALSE,
             size = 2, shape = 21, fill = "black", color = "black", stroke = 0.5) +
  geom_text(data = letters_cl, aes(x = Climate, y = y_pos, label = group_letters),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  facet_grid(Response ~ ., scales = "free_y") +
  scale_fill_manual(values = cc) +
  labs(x = "Climatic conditions", y = "Observed value") +
  theme_facet +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold", size = 10, color = "black"))

## 6.3 Species × Climate violins (interaction)
df_violin_SC <- dplyr::bind_rows(
  data_biomass |>
    dplyr::select(Mesocosm, Species, Climate, Value = Species_biomass) |>
    dplyr::mutate(Response = "Species biomass"),
  data_biomass |>
    dplyr::select(Mesocosm, Species, Climate, Value = Proportion) |>
    dplyr::mutate(Response = "Biomass proportion"),
  data_traits |>
    dplyr::select(Mesocosm, Species, Climate, Value = SLA) |>
    dplyr::mutate(Response = "SLA"),
  data_traits |>
    dplyr::select(Mesocosm, Species, Climate, Value = LDMC) |>
    dplyr::mutate(Response = "LDMC")
) |>
  dplyr::filter(!is.na(Value)) |>
  dplyr::mutate(Species  = factor(Species, levels = species_levels),
                Climate  = factor(Climate, levels = c("Current", "Future")),
                Response = factor(Response,
                                  levels = c("Species biomass", "Biomass proportion",
                                             "SLA", "LDMC")))

df_violin_SC <- df_violin_SC |>
  dplyr::mutate(Climate_num = as.numeric(Climate))

iqr_SC <- df_violin_SC |>
  dplyr::group_by(Response, Climate, Climate_num, Species) |>
  dplyr::summarise(median = median(Value, na.rm = TRUE),
                   q1     = quantile(Value, 0.25, na.rm = TRUE),
                   q3     = quantile(Value, 0.75, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::mutate(dodge_width   = 0.8,
                offset        = ifelse(Species == "Festuca", -0.5, 0.5) * (dodge_width / 2),
                Climate_dodge = Climate_num + offset)

letters_SC <- df_violin_SC |>
  dplyr::group_by(Response) |>
  dplyr::group_modify(~ {
    resp <- .y$Response
    if (resp %in% c("SLA", "LDMC")) {
      m <- lme4::lmer(Value ~ Species * Climate + (1 | Mesocosm), data = .x)
    } else {
      m <- stats::lm(Value ~ Species * Climate, data = .x)
    }
    emm <- emmeans::emmeans(m, ~ Species:Climate)
    cl  <- multcomp::cld(emm, adjust = "sidak", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Species, Climate, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(Climate_num = as.numeric(Climate),
                dodge_width = 0.8,
                offset      = ifelse(Species == "Festuca", -0.5, 0.5) * (dodge_width / 2),
                Climate_dodge = Climate_num + offset) |>
  dplyr::left_join(df_violin_SC |>
                     dplyr::group_by(Response) |>
                     dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop"),
                   by = "Response") |>
  dplyr::mutate(y_pos = y_max * 1.10)

p_violin_SC <- ggplot(df_violin_SC,
                      aes(x = Climate, y = Value, fill = Species)) +
  geom_violin(position = position_dodge(width = 0.8),
              trim = FALSE, alpha = 0.7, color = NA) +
  geom_point(data = iqr_SC, aes(x = Climate, y = median, group = Species),
             position = position_dodge(width = 0.8), inherit.aes = FALSE,
             size = 2, shape = 21, color = "black", fill = "black", stroke = 0.5) +
  geom_segment(data = iqr_SC, aes(x = Climate_dodge, xend = Climate_dodge, y = q1, yend = q3),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_text(data = letters_SC, aes(x = Climate_dodge, y = y_pos, label = group_letters),
            inherit.aes = FALSE, fontface = "bold", size = 3) +
  facet_grid(Response ~ ., scales = "free_y") +
  scale_fill_manual(values = cs) +
  labs(x = "Climatic conditions", y = "Observed value", fill = "Plant species") +
  theme_facet +
  theme(legend.position = "top",
        legend.text = element_text(face = "bold.italic", size = 10),
        axis.text.x  = element_text(face = "bold", size = 11))


# -------------------------------------------------------------------
# 7) Final combined figures & saving
# -------------------------------------------------------------------

# 7.1 Species-focused figure (Species violins + Moisture × Species & Climate)
legend_MC <- ggpubr::as_ggplot(ggpubr::get_legend(p_raw_sig_clim))
p_MC_noleg <- p_raw_sig_clim + theme(legend.position = "none")

f_species_inner <- ggpubr::ggarrange(final_plot_MP + theme(legend.position = "none"),
                                     p_MC_noleg, nrow = 1, widths = c(1, 1.2),
                                     labels = c("(b)", "(c)"), font.label = list(size = 11))

f_species <- ggpubr::ggarrange(p_violin_sp, f_species_inner, legend_MC, nrow = 3,
                               heights = c(1, 1.8, 0.12), labels = c("(a)", "", ""),
                               font.label = list(size = 11))

ggsave(file.path(fig_dir, "Plant_performance_Species.png"), ## Figure 1
       f_species, bg = "white", width = 14, height = 8, units = "in", dpi = 600)

# 7.2 Moisture × Climate figure (M-only + M×C)
f_MxC <- ggpubr::ggarrange(final_plot_M, final_plot_MC, nrow = 1,
                           widths = c(1, 1), labels = c("(a)", "(b)"),
                           font.label = list(size = 11))

ggsave(file.path(fig_dir, "Plant_performance_MxC.png"), ## Figure S6
       f_MxC, bg = "white", width = 12, height = 8, units = "in", dpi = 600)

# 7.3 Climate-focused figure (Climate-only + Species×Climate violins)
f_climate <- ggpubr::ggarrange(p_violin_cl, p_violin_SC, nrow = 1,
                               widths = c(1, 1.8), labels = c("(a)", "(b)"),
                               font.label = list(size = 11))

ggsave(file.path(fig_dir, "Plant_performance_Climate.png"), ## Figure S7
       f_climate, bg = "white", width = 12, height = 8, units = "in", dpi = 600)


