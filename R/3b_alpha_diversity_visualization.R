##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/3b_alpha_diversity_visualization.R
##
## Tasks:
##  1) Visualise microbial ASV richness vs Moisture:
##       - Species × Moisture (Bacteria / Fungi)
##       - Climate × Moisture (Bacteria / Fungi)
##       - Overall Moisture trends across compartments
##       - Species × Moisture × Climate (Bacteria / Fungi)
##  2) Violin plots:
##       - Species × Climate (by compartment)
##       - Species-only
##       - Climate-only
##  3) Export composite figures:
##       - Microbial Richness_Species.png
##       - Microbial Richness_M+C.png
##       - Microbial Richness_Climate.png
##
## Notes:
##  - Uses project-wide utilities in R/utils_packages.R and R/utils_theme.R
##  - Uses add_drought_zones() and project_theme_facet() for consistent styling

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
data_dir    <- file.path(root, "Data", "3_Alpha diversity")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

theme_facet <- project_theme_facet()
cs <- species_cols   # from utils_theme.R (Festuca/Lolium)
cc <- climate_cols   # from utils_theme.R (Current/Future)

# -------------------------------------------------------------------
# 1) Import and reshape richness data
# -------------------------------------------------------------------
rich <- read.csv(file.path(data_dir, "Richness_all.csv"),
                 check.names = FALSE)

# Factor levels consistent with biomass / other scripts
species_levels <- levels(factor(rich$Species))
climate_levels <- levels(factor(rich$Climate))

rich <- rich |>
  dplyr::mutate(Species = factor(Species, levels = species_levels),
                Climate = factor(Climate, levels = climate_levels))

# Long format for richness (Domain × Compartment)
rich_long <- rich |>
  tidyr::pivot_longer(cols = dplyr::starts_with("Richness_"),
                      names_to = c("Domain", "Compartment"),
                      names_pattern = "Richness_(Bac|Fun)_(Leaf|Root|Rhizo)",
                      values_to = "Richness") |>
  dplyr::mutate(Domain = dplyr::recode(Domain, "Bac" = "Bacteria", "Fun" = "Fungi"),
                Compartment = dplyr::recode(Compartment,
                                            "Leaf"  = "Leaf endosphere",
                                            "Root"  = "Root endosphere",
                                            "Rhizo" = "Rhizosphere"),
                Compartment = factor(Compartment, levels = c("Leaf endosphere",
                                                             "Root endosphere",
                                                             "Rhizosphere")))

# -------------------------------------------------------------------
# 2) Helper functions
# -------------------------------------------------------------------

# Drought labels strip (top panel) using same colors as add_drought_zones()
make_drought_top_strip <- function(x_min, x_max, thr1 = 400, thr2 = 1000,
                                   labels = c("Severe drought",
                                              "Moderate drought",
                                              "Well-watered")) {
  mid_severe   <- (x_min + thr1) / 2
  mid_moderate <- (thr1 + thr2) / 2
  mid_wet      <- (thr2 + x_max) / 2
  
  lab_df <- tibble::tibble(x     = c(mid_severe, mid_moderate, mid_wet),
                           y     = 0.5, label = labels)
  
  ggplot2::ggplot(lab_df, ggplot2::aes(x = x, y = y, label = label)) +
    ggplot2::annotate("rect", xmin = x_min, xmax = thr1, ymin = 0, ymax = 1,
                      fill = "#f9d7c2", alpha = 0.4) +
    ggplot2::annotate("rect", xmin = thr1, xmax = thr2, ymin = 0, ymax = 1,
                      fill = "#fffec0", alpha = 0.4) +
    ggplot2::annotate("rect", xmin = thr2, xmax = x_max, ymin = 0, ymax = 1,
                      fill = "#d9eaf7", alpha = 0.4) +
    ggplot2::geom_text(fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 1, r = 15, b = -10, l = 15))
}

# -------------------------------------------------------------------
# 3) Species × Moisture (by domain, 3 compartments)
# -------------------------------------------------------------------

# Helper: p-value → stars (from utils_theme.R; here we just call p_to_stars())

make_domain_plot_species <- function(df_domain, y_lab = "ASV richness") {
  
  # Slopes per Compartment × Species
  slopes <- df_domain |>
    dplyr::group_by(Compartment, Species) |>
    dplyr::group_modify(~ {
      m  <- stats::lm(Richness ~ Moisture, data = .x)
      sm <- broom::tidy(m)
      sm |>
        dplyr::filter(term == "Moisture") |>
        dplyr::select(estimate, p.value)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                  stars     = p_to_stars(p.value),
                  label     = paste0("Slope = ",
                                     formatC(estimate, format = "f", digits = 3), " ",
                                     stars))
  
  df2 <- df_domain |>
    dplyr::left_join(slopes, by = c("Compartment", "Species"))
  
  x_min <- min(df2$Moisture, na.rm = TRUE)
  x_max <- max(df2$Moisture, na.rm = TRUE)
  thr1  <- 400
  thr2  <- 1000
  
  panel_range <- df2 |>
    dplyr::group_by(Compartment) |>
    dplyr::summarise(x_min = min(Moisture, na.rm = TRUE),
                     y_max = max(Richness, na.rm = TRUE),
                     .groups = "drop")
  
  annot_slopes <- slopes |>
    dplyr::left_join(panel_range, by = "Compartment") |>
    dplyr::group_by(Compartment) |>
    dplyr::arrange(Species) |>
    dplyr::mutate(x = x_max * 0.64,
                  y = y_max * (0.95 - 0.10 * (dplyr::row_number() - 1))) |>
    dplyr::ungroup()
  
  p_raw <- ggplot2::ggplot(df2, ggplot2::aes(x = Moisture, y = Richness, color = Species)) +
    add_drought_zones(x_min, x_max, thr1, thr2) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE,
                         ggplot2::aes(group = Species, fill = Species, linetype = line_type),
                         alpha = 0.15, color = NA, linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE,
                         ggplot2::aes(group = Species, color = Species, linetype = line_type),
                         linewidth = 1.0) +
    ggplot2::facet_grid(Compartment ~ ., scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = cs) +
    ggplot2::scale_fill_manual(values = cs) +
    ggplot2::scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed"), drop = FALSE) +
    ggplot2::guides(linetype = ggplot2::guide_legend(title = "Significance",
                                                     override.aes = list(color = "black",
                                                                         linewidth = 0.5))) +
    ggplot2::labs(x = "Moisture (watering amount, ml)", y = y_lab,
                  color   = "Plant species", linetype = "Significance") +
    theme_facet +
    ggplot2::geom_text(data = annot_slopes,
                       ggplot2::aes(x = x, y = y, label = label, color = Species),
                       inherit.aes = FALSE, hjust = 0, vjust = 1,
                       size = 3.5, fontface = "bold", show.legend = FALSE)
  
  p_top <- make_drought_top_strip(x_min, x_max, thr1, thr2)
  
  p_top / p_raw + patchwork::plot_layout(heights = c(0.06, 1))
}

rich_bac <- rich_long |>
  dplyr::filter(Domain == "Bacteria")
rich_fun <- rich_long |>
  dplyr::filter(Domain == "Fungi")

plot_bacteria <- make_domain_plot_species(rich_bac,
                                          y_lab = "Bacterial ASV richness")
plot_fungi    <- make_domain_plot_species(rich_fun,
                                          y_lab = "Fungal ASV richness")

plot_bacteria_noleg <- plot_bacteria + ggplot2::theme(legend.position = "none")
plot_fungi_noleg    <- plot_fungi    + ggplot2::theme(legend.position = "none")

plot_bac_top <- plot_bacteria +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(color   = ggplot2::guide_legend(nrow = 1, order = 1),
                  linetype = ggplot2::guide_legend(
                    nrow = 1, order = 2,
                    override.aes = list(color = "black", linewidth = 0.5)))

legend_bac <- ggplotify::as.ggplot(ggpubr::get_legend(plot_bac_top))

p3 <- ggpubr::ggarrange(plot_bacteria_noleg, plot_fungi_noleg, ncol = 1,
                        heights = c(1, 1.4), labels = c("(b)", "(c)"),
                        font.label = list(size = 11))

# -------------------------------------------------------------------
# 4) Climate × Moisture (by domain, 3 compartments)
# -------------------------------------------------------------------

make_domain_plot_climate <- function(df_domain, y_lab = "ASV richness") {
  
  slopes <- df_domain |>
    dplyr::group_by(Compartment, Climate) |>
    dplyr::group_modify(~ {
      m  <- stats::lm(Richness ~ Moisture, data = .x)
      sm <- broom::tidy(m)
      sm |>
        dplyr::filter(term == "Moisture") |>
        dplyr::select(estimate, p.value)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                  stars     = p_to_stars(p.value),
                  label     = paste0("Slope = ",
                                     formatC(estimate, format = "f", digits = 3), " ",
                                     stars))
  
  df2 <- df_domain |>
    dplyr::left_join(slopes, by = c("Compartment", "Climate"))
  
  x_min <- min(df2$Moisture, na.rm = TRUE)
  x_max <- max(df2$Moisture, na.rm = TRUE)
  thr1  <- 400
  thr2  <- 1000
  
  panel_range <- df2 |>
    dplyr::group_by(Compartment) |>
    dplyr::summarise(x_min = min(Moisture, na.rm = TRUE),
                     y_max = max(Richness, na.rm = TRUE), .groups = "drop")
  
  annot_slopes <- slopes |>
    dplyr::left_join(panel_range, by = "Compartment") |>
    dplyr::group_by(Compartment) |>
    dplyr::arrange(Climate) |>
    dplyr::mutate(x = x_max * 0.64,
                  y = y_max * (0.95 - 0.10 * (dplyr::row_number() - 1))) |>
    dplyr::ungroup()
  
  p_raw <- ggplot2::ggplot(df2,
                           ggplot2::aes(x = Moisture,
                                        y = Richness,
                                        color = Climate)) +
    add_drought_zones(x_min, x_max, thr1, thr2) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE,
                         ggplot2::aes(group = Climate, fill = Climate, linetype = line_type),
                         alpha = 0.15, color = NA, linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE,
                         ggplot2::aes(group = Climate, color = Climate, linetype = line_type),
                         linewidth = 1.0) +
    ggplot2::facet_grid(Compartment ~ ., scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed"), drop = FALSE) +
    ggplot2::guides(linetype = ggplot2::guide_legend(title = "Significance",
                                                     override.aes = list(color = "black", linewidth = 0.5))) +
    ggplot2::labs(x = "Moisture (watering amount, ml)", y = y_lab,
                  color   = "Climatic conditions", linetype = "Significance") +
    theme_facet +
    ggplot2::geom_text(data = annot_slopes,
                       ggplot2::aes(x = x, y = y, label = label, color = Climate),
                       inherit.aes = FALSE, hjust = 0, vjust = 1,
                       size = 3.5, fontface = "bold", show.legend = FALSE)
  
  p_top <- make_drought_top_strip(x_min, x_max, thr1, thr2)
  
  p_top / p_raw + patchwork::plot_layout(heights = c(0.06, 1))
}

plot_bacteria_clim <- make_domain_plot_climate(rich_bac,
                                               y_lab = "Bacterial ASV richness")
plot_fungi_clim    <- make_domain_plot_climate(rich_fun,
                                               y_lab = "Fungal ASV richness")

plot_bacteria_clim_noleg <- plot_bacteria_clim + ggplot2::theme(legend.position = "none")
plot_fungi_clim_noleg    <- plot_fungi_clim    + ggplot2::theme(legend.position = "none")

plot_bac_clim_top <- plot_bacteria_clim +
  ggplot2::theme(legend.position = "top") +
  ggplot2::guides(color   = ggplot2::guide_legend(nrow = 1, order = 1),
                  linetype = ggplot2::guide_legend(
                    nrow = 1, order = 2,
                    override.aes = list(color = "black", linewidth = 0.5)))

legend_bac_clim <- ggplotify::as.ggplot(ggpubr::get_legend(plot_bac_clim_top))

pc3 <- ggpubr::ggarrange(plot_bacteria_clim_noleg,
                         plot_fungi_clim_noleg,
                         legend_bac_clim,
                         ncol = 1, labels = c("(b)", "(c)", " "),
                         heights = c(1, 1.4, 0.05),
                         font.label = list(size = 11))

# -------------------------------------------------------------------
# 5) Overall Moisture (pooled within panel = Domain × Compartment)
# -------------------------------------------------------------------
rich_both <- rich_long |>
  dplyr::filter(!(Domain == "Bacteria" & Compartment == "Leaf endosphere")) |>
  dplyr::mutate(Panel = interaction(Domain, Compartment, sep = " – ")) |>
  dplyr::mutate(Panel = dplyr::recode(Panel,
                                      "Fungi – Leaf endosphere"    = "Fungi – Leaf",
                                      "Fungi – Root endosphere"    = "Fungi – Root",
                                      "Fungi – Rhizosphere"        = "Fungi – Rhizo",
                                      "Bacteria – Root endosphere" = "Bacteria – Root",
                                      "Bacteria – Rhizosphere"     = "Bacteria – Rhizo"))

panel_levels <- c("Bacteria – Root", "Bacteria – Rhizo",
                  "Fungi – Leaf", "Fungi – Root", "Fungi – Rhizo")

rich_both <- rich_both |>
  dplyr::mutate(Panel = factor(Panel, levels = panel_levels))

slopes_both <- rich_both |>
  dplyr::group_by(Panel) |>
  dplyr::group_modify(~ {
    m  <- stats::lm(Richness ~ Moisture, data = .x)
    sm <- broom::tidy(m)
    sm |>
      dplyr::filter(term == "Moisture") |>
      dplyr::select(estimate, p.value)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                stars     = p_to_stars(p.value),
                label     = paste0("Slope = ",
                                   formatC(estimate, format = "f", digits = 3), " ",
                                   stars))

rich_both2 <- rich_both |>
  dplyr::left_join(slopes_both, by = "Panel")

x_min <- min(rich_both2$Moisture, na.rm = TRUE)
x_max <- max(rich_both2$Moisture, na.rm = TRUE)
thr1  <- 400
thr2  <- 1000

panel_range <- rich_both2 |>
  dplyr::group_by(Panel) |>
  dplyr::summarise(x_min = min(Moisture, na.rm = TRUE),
                   y_max = max(Richness, na.rm = TRUE), .groups = "drop")

annot_slopes_both <- slopes_both |>
  dplyr::left_join(panel_range, by = "Panel") |>
  dplyr::group_by(Panel) |>
  dplyr::mutate(x = x_max * 0.64, y = y_max * 0.95) |>
  dplyr::ungroup()

p_overall_both <- ggplot2::ggplot(rich_both2,
                                  ggplot2::aes(x = Moisture, y = Richness, linetype = line_type)) +
  add_drought_zones(x_min, x_max, thr1, thr2) +
  ggplot2::geom_point(alpha = 0.4, size = 1.5, color = "grey60") +
  ggplot2::geom_smooth(method = "lm", se = TRUE,
                       ggplot2::aes(linetype = line_type, group = Panel),
                       linewidth = 0.4, color = NA, fill = "grey50",
                       alpha = 0.2, show.legend = FALSE) +
  ggplot2::geom_smooth(method = "lm", se = FALSE,
                       ggplot2::aes(linetype = line_type, group = Panel),
                       linewidth = 1.0, color = "blue", show.legend = TRUE) +
  ggplot2::facet_grid(Panel ~ ., scales = "free_y") +
  ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  ggplot2::scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed"),
                                 name   = "Significance") +
  ggplot2::guides(linetype = ggplot2::guide_legend(
    title = "Significance", nrow = 1, override.aes = list(color = "black", linewidth = 0.5))) +
  ggplot2::labs(x = "Moisture (watering amount, ml)", y = "ASV richness") +
  theme_facet +
  ggplot2::theme(legend.position        = c(0.35, 0.03),
                 legend.title.position  = "left",
                 legend.background      = ggplot2::element_rect(
                   fill = "white", color = "black", linewidth = 0.3, linetype = "dotted")) +
  ggplot2::geom_text(data = annot_slopes_both, ggplot2::aes(x = x, y = y, label = label),
                     inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5, fontface = "bold")

p_top_both <- make_drought_top_strip(x_min, x_max, thr1, thr2)

overall_bac_fun_5 <- p_top_both / p_overall_both +
  patchwork::plot_layout(heights = c(0.06, 1))

# -------------------------------------------------------------------
# 6) Species × Moisture × Climate (by domain)
# -------------------------------------------------------------------

make_domain_plot_species_climate <- function(df_domain,
                                             y_lab = "ASV richness") {
  
  slopes <- df_domain |>
    dplyr::group_by(Compartment, Species, Climate) |>
    dplyr::group_modify(~ {
      m  <- stats::lm(Richness ~ Moisture, data = .x)
      sm <- broom::tidy(m)
      sm |>
        dplyr::filter(term == "Moisture") |>
        dplyr::select(estimate, p.value)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(line_type = ifelse(p.value < 0.05, "sig", "nonsig"),
                  stars     = p_to_stars(p.value),
                  label     = paste0("Slope = ",
                                     formatC(estimate, format = "f", digits = 3), " ",
                                     stars))
  
  df2 <- df_domain |>
    dplyr::left_join(slopes, by = c("Compartment", "Species", "Climate"))
  
  x_min <- min(df2$Moisture, na.rm = TRUE)
  x_max <- max(df2$Moisture, na.rm = TRUE)
  thr1  <- 400
  thr2  <- 1000
  
  panel_range <- df2 |>
    dplyr::group_by(Compartment, Climate) |>
    dplyr::summarise(x_min = min(Moisture, na.rm = TRUE),
                     y_max = max(Richness, na.rm = TRUE), .groups = "drop")
  
  annot_slopes <- slopes |>
    dplyr::left_join(panel_range, by = c("Compartment", "Climate")) |>
    dplyr::group_by(Compartment, Climate) |>
    dplyr::arrange(Species) |>
    dplyr::mutate(x = x_max * 0.64, y = y_max * (0.95 - 0.10 * (dplyr::row_number() - 1))) |>
    dplyr::ungroup()
  
  p_raw <- ggplot2::ggplot(df2, ggplot2::aes(x = Moisture, y = Richness, color = Species)) +
    add_drought_zones(x_min, x_max, thr1, thr2) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE,
                         ggplot2::aes(group   = interaction(Species, Climate),
                                      fill    = Species, linetype = line_type),
                         alpha = 0.15, color = NA, linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE,
                         ggplot2::aes(group   = interaction(Species, Climate),
                                      color  = Species, linetype = line_type),
                         linewidth = 1.0) +
    ggplot2::facet_grid(Compartment ~ Climate, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = cs) +
    ggplot2::scale_fill_manual(values = cs) +
    ggplot2::scale_linetype_manual(values = c(sig = "solid", nonsig = "dashed"), drop = FALSE) +
    ggplot2::guides(linetype = ggplot2::guide_legend(
      title = "Significance", override.aes = list(color = "black", linewidth = 0.5))) +
    ggplot2::labs(x = "Moisture (watering amount, ml)", y = y_lab,
                  color   = "Plant species", linetype = "Significance") +
    theme_facet +
    ggplot2::geom_text(data = annot_slopes,
                       ggplot2::aes(x = x, y = y, label = label, color = Species),
                       inherit.aes = FALSE, hjust = 0, vjust = 1,
                       size = 3.5, fontface = "bold", show.legend = FALSE)
  p_raw
}

plot_bacteria_SC <- make_domain_plot_species_climate(rich_bac,
                                                     y_lab = "Bacterial ASV richness")

plot_fungi_SC <- make_domain_plot_species_climate(rich_fun,
                                                  y_lab = "Fungal ASV richness")

plot_b_noleg <- plot_bacteria_SC + ggplot2::theme(legend.position = "none")
plot_f_noleg <- plot_fungi_SC    + ggplot2::theme(legend.position = "none")

psc4 <- ggpubr::ggarrange(plot_b_noleg, plot_f_noleg, ncol = 1, heights = c(1, 1.4),
                          labels = c("(d)", "(e)"), font.label = list(size = 11))

psc5 <- ggpubr::ggarrange(p3, psc4, nrow = 1, widths = c(1, 1.3))
psc6 <- ggpubr::ggarrange(psc5, legend_bac, ncol = 1, heights = c(1, 0.05))

# -------------------------------------------------------------------
# 7) Species × Climate violins (5 panels, Species × Climate)
# -------------------------------------------------------------------
rich_both_sc <- rich_long |>
  dplyr::filter(!(Domain == "Bacteria" & Compartment == "Leaf endosphere")) |>
  dplyr::mutate(Panel = interaction(Domain, Compartment, sep = " – ")) |>
  dplyr::mutate(Panel = dplyr::recode(Panel,
                                      "Bacteria – Root endosphere" = "Bacteria – Root",
                                      "Bacteria – Rhizosphere"     = "Bacteria – Rhizo",
                                      "Fungi – Leaf endosphere"    = "Fungi – Leaf",
                                      "Fungi – Root endosphere"    = "Fungi – Root",
                                      "Fungi – Rhizosphere"        = "Fungi – Rhizo"),
                Panel   = factor(Panel, levels = panel_levels),
                Species = factor(Species), Climate = factor(Climate))

df_sc <- rich_both_sc |>
  dplyr::filter(!is.na(Richness)) |>
  dplyr::rename(Value = Richness) |>
  dplyr::mutate(Climate_num = as.numeric(Climate))

iqr_sc <- df_sc |>
  dplyr::group_by(Panel, Climate, Climate_num, Species) |>
  dplyr::summarise(median = stats::median(Value, na.rm = TRUE),
                   q1     = stats::quantile(Value, 0.25, na.rm = TRUE),
                   q3     = stats::quantile(Value, 0.75, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(dodge_width   = 0.8,
                offset        = ifelse(Species == "Festuca", -0.5, 0.5) * (dodge_width / 2),
                Climate_dodged = Climate_num + offset)

letters_sc <- df_sc |>
  dplyr::group_by(Panel) |>
  dplyr::group_modify(~ {
    m   <- lme4::lmer(Value ~ Species * Climate + (1 | Mesocosm), data = .x)
    emm <- emmeans::emmeans(m, ~ Species:Climate)
    cl  <- multcomp::cld(emm, adjust = "sidak", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Species, Climate, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    df_sc |>
      dplyr::group_by(Panel) |>
      dplyr::summarise(
        y_pos = max(Value, na.rm = TRUE) * 1.05,
        .groups = "drop"
      ),
    by = "Panel"
  )

pos_dodge <- ggplot2::position_dodge(width = 0.8)

p_sc <- ggplot2::ggplot(df_sc, ggplot2::aes(x = Climate, y = Value, fill = Species)) +
  ggplot2::geom_violin(position = pos_dodge, trim = FALSE, alpha = 0.7, color = NA) +
  ggplot2::geom_point(data = iqr_sc, ggplot2::aes(x = Climate, y = median, group = Species),
                      position = pos_dodge, inherit.aes = FALSE, size = 2, shape = 21,
                      fill = "black", color = "black", stroke = 0.5) +
  ggplot2::geom_segment(data = iqr_sc, ggplot2::aes(
    x = Climate_dodged, xend = Climate_dodged, y = q1, yend = q3),
    inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  ggplot2::facet_grid(Panel ~ ., scales = "free_y") +
  ggplot2::scale_fill_manual(values = cs) +
  ggplot2::labs(x = "Climatic conditions", y = "ASV richness", fill = "Plant species") +
  ggplot2::geom_text(data = letters_sc,
                     ggplot2::aes(x = Climate, y = y_pos, label = group_letters, group = Species),
                     position = pos_dodge, inherit.aes = FALSE,
                     fontface = "bold", size = 3) +
  theme_facet +
  ggplot2::theme(legend.position = "top",
                 axis.text.x     = ggplot2::element_text(face = "bold", size = 10))

# -------------------------------------------------------------------
# 8) Species-only violins
# -------------------------------------------------------------------
df_sp <- rich_both_sc |>
  dplyr::filter(!is.na(Richness)) |>
  dplyr::rename(Value = Richness)

iqr_sp <- df_sp |>
  dplyr::group_by(Panel, Species) |>
  dplyr::summarise(median = stats::median(Value, na.rm = TRUE),
                   q1     = stats::quantile(Value, 0.25, na.rm = TRUE),
                   q3     = stats::quantile(Value, 0.75, na.rm = TRUE),
                   .groups = "drop")

letters_sp <- df_sp |>
  dplyr::group_by(Panel) |>
  dplyr::group_modify(~ {
    m   <- lme4::lmer(Value ~ Species + (1 | Mesocosm), data = .x)
    emm <- emmeans::emmeans(m, ~ Species)
    cl  <- multcomp::cld(emm, adjust = "none", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Species, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    df_sp |>
      dplyr::group_by(Panel) |>
      dplyr::summarise(y_pos = max(Value, na.rm = TRUE) * 1.05, .groups = "drop"),
    by = "Panel")

p_species <- ggplot2::ggplot(df_sp, ggplot2::aes(x = Species, y = Value, fill = Species)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  ggplot2::geom_point(data = iqr_sp, ggplot2::aes(x = Species, y = median),
                      inherit.aes = FALSE, size = 2, shape = 21,
                      fill = "black", color = "black", stroke = 0.5) +
  ggplot2::geom_segment(data = iqr_sp,
                        ggplot2::aes(x = Species, xend = Species, y = q1, yend = q3),
                        inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  ggplot2::facet_grid(Panel ~ ., scales = "free_y") +
  ggplot2::scale_fill_manual(values = cs) +
  ggplot2::labs(x = "Plant species", y = "ASV richness") +
  ggplot2::geom_text(data = letters_sp,
                     ggplot2::aes(x = Species, y = y_pos, label = group_letters),
                     inherit.aes = FALSE, fontface = "bold", size = 3) +
  theme_facet +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(face  = "bold.italic",
                                                     size  = 10,
                                                     color = "black"))

# -------------------------------------------------------------------
# 9) Climate-only violins
# -------------------------------------------------------------------
df_cl <- rich_both_sc |>
  dplyr::filter(!is.na(Richness)) |>
  dplyr::rename(Value = Richness)

iqr_cl <- df_cl |>
  dplyr::group_by(Panel, Climate) |>
  dplyr::summarise(median = stats::median(Value, na.rm = TRUE),
                   q1     = stats::quantile(Value, 0.25, na.rm = TRUE),
                   q3     = stats::quantile(Value, 0.75, na.rm = TRUE),
                   .groups = "drop")

letters_cl <- df_cl |>
  dplyr::group_by(Panel) |>
  dplyr::group_modify(~ {
    m   <- lme4::lmer(Value ~ Climate + (1 | Mesocosm), data = .x)
    emm <- emmeans::emmeans(m, ~ Climate)
    cl  <- multcomp::cld(emm, adjust = "none", Letters = letters) |>
      as.data.frame()
    cl |>
      dplyr::select(Climate, .group) |>
      dplyr::rename(group_letters = .group)
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    df_cl |>
      dplyr::group_by(Panel) |>
      dplyr::summarise(
        y_pos = max(Value, na.rm = TRUE) * 1.05,
        .groups = "drop"
      ),
    by = "Panel"
  )

p_climate <- ggplot2::ggplot(df_cl, ggplot2::aes(x = Climate, y = Value, fill = Climate)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  ggplot2::geom_point(data = iqr_cl,
                      ggplot2::aes(x = Climate, y = median),
                      inherit.aes = FALSE, size = 2, shape = 21,
                      fill = "black", color = "black", stroke = 0.5) +
  ggplot2::geom_segment(data = iqr_cl,
                        ggplot2::aes(x = Climate, xend = Climate, y = q1, yend = q3),
                        inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  ggplot2::facet_grid(Panel ~ ., scales = "free_y") +
  ggplot2::scale_fill_manual(values = cc) +
  ggplot2::labs(x    = "Climatic conditions", y    = "ASV richness",
                fill = "Climatic conditions") +
  ggplot2::geom_text(data = letters_cl,
                     ggplot2::aes(x = Climate, y = y_pos, label = group_letters),
                     inherit.aes = FALSE, fontface = "bold", size = 3) +
  theme_facet +
  ggplot2::theme(legend.position = "none",
                 axis.text.x     = ggplot2::element_text(face  = "bold",
                                                         size  = 10,
                                                         color = "black"))

# -------------------------------------------------------------------
# 10) Final composite figures & saving
# -------------------------------------------------------------------

# 10.1 Species-focused figure (violins + species/climate interactions)
f4 <- ggpubr::ggarrange(p_species, psc6, nrow = 1, widths = c(1, 5),
                        labels = c("(a)", " "), font.label = list(size = 11),
                        label.x = -0.01)

ggplot2::ggsave(file.path(results_dir, "Microbial Richness_Species.png"), ## Figure 2
                f4, bg = "white", width = 16, height = 10, units = "in", dpi = 600)

# 10.2 Moisture + Climate figure (overall + climate × moisture)
f5 <- ggpubr::ggarrange(overall_bac_fun_5, pc3, nrow = 1, widths = c(1, 1),
                        labels = c("(a)", " "), font.label = list(size = 11),
                        label.x = -0.01)

ggplot2::ggsave(file.path(results_dir, "Microbial Richness_M+C.png"), ## Figure S8
                f5, bg = "white", width = 12, height = 9, units = "in", dpi = 600)

# 10.3 Climate-focused figure (climate-only & Species × Climate violins)
f6 <- ggpubr::ggarrange(p_climate, p_sc, nrow = 1, widths = c(1, 1.8),
                        labels = c("(a)", "(b)"), font.label = list(size = 11),
                        label.x = -0.01)

ggplot2::ggsave(file.path(results_dir, "Microbial Richness_Climate.png"), ## Figure S9
                f6, bg = "white", width = 8, height = 8, units = "in", dpi = 600)

