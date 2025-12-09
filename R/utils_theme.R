##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## R/utils_theme.R

project_theme_classic <- function() {
  theme_classic() +
    theme(axis.text   = element_text(face = "bold", size = 10),
          axis.title  = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 12),
          legend.text  = element_text(face = "bold", size = 12),
          legend.position = "top")
}

project_theme_facet <- function() {
  theme_classic() +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 10),
          legend.text  = element_text(face = "bold", size = 10),
          strip.text   = element_text(face = "bold", size = 11),
          strip.background = element_rect(color = "black", fill = "gray98", linewidth = 0.5),
          axis.text   = element_text(size = 9),
          axis.title  = element_text(face = "bold", size = 11),
          panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

# Color palettes
species_cols <- c(Festuca = "#E59400", Lolium = "#198C19")
climate_cols <- c(Current = "#377EB8", Future = "#E94849")

# Simple helper to convert p-values to stars
p_to_stars <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***",
                   p < 0.01  ~ "**",
                   p < 0.05  ~ "*",
                   TRUE      ~ "ns")
}

# Drought-zone background rects (for Moisture plots)
add_drought_zones <- function(x_min, x_max,
                              thr1 = 400, thr2 = 1000,
                              severe_col   = "#f9d7c2",
                              moderate_col = "#fffec0",
                              wet_col      = "#d9eaf7") {
  list(annotate("rect", xmin = x_min, xmax = thr1, ymin = -Inf, ymax = Inf,
             fill = severe_col, alpha = 0.4),
       annotate("rect", xmin = thr1, xmax = thr2, ymin = -Inf, ymax = Inf,
             fill = moderate_col, alpha = 0.4),
       annotate("rect", xmin = thr2, xmax = x_max, ymin = -Inf, ymax = Inf,
             fill = wet_col, alpha = 0.4))
}
