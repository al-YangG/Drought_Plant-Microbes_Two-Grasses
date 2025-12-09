##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/5_variance_explained.R
##
## NOTE:
##   The file "Variance Explained_Plant+Microbes.csv" is a manually
##   prepared combined dataset that merges variance-explained outputs from:
##     (1) Plant performance models
##     (2) Microbial richness models (alpha diversity)
##     (3) Microbial community composition models (RDA)
##
##   This script ONLY reads the combined table and produces a unified
##   heatmap of variance explained (%) across all response variables.
##
## Tasks:
##   1) Read combined variance explained table
##   2) Harmonize factor/response names
##   3) Plot heatmap of variance explained (%) across groups
##   4) Save as results/Variance_explained_heatmap.png

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
data_dir    <- file.path(root, "Data", "5_Variation explained")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Import combined variance-explained table
# -------------------------------------------------------------------
df <- read.csv(file.path(data_dir, "Variance Explained_Plant+Microbes.csv"),
               na.strings = c("", "NA"), check.names = FALSE)

# -------------------------------------------------------------------
# 2) Recode labels for compact plotting
# -------------------------------------------------------------------
df <- df |>
  dplyr::mutate(Group_full = Group,
                Group_short = dplyr::recode(Group,
                                            "Plant performance"               = "Plant",
                                            "Microbial richness"              = "M_Rich",
                                            "Microbial community composition" = "M_Comp",
                                            .default = Group),
    Response_short = Response |>
      stringr::str_replace_all("Bacteria", "Bac") |>
      stringr::str_replace_all("Fungi", "Fun") |>
      stringr::str_replace_all("Rhizosphere", "Rhizo") |>
      stringr::str_replace_all("Root endosphere", "Root") |>
      stringr::str_replace_all("Leaf endosphere", "Leaf") |>
      stringr::str_replace_all("Species biomass", "Biomass"),
    Row_short = paste0(Group_short, " | ", Response_short))

# -------------------------------------------------------------------
# 3) Ordering for axes and facets
# -------------------------------------------------------------------
factor_levels <- c("Moisture (M)", "Climatic conditions (C)", "Plant species (P)",
                   "M:C", "M:P", "C:P", "M:C:P")

factor_labels <- c("Moisture (M)", "Climatic conditions (C)", "Plant species (P)",
                   "M × C", "M × P", "C × P", "M × C × P")

row_levels <- c("Plant | Biomass", "Plant | SLA", "Plant | LDMC",
                "M_Rich | Bac_Root", "M_Rich | Bac_Rhizo",
                "M_Rich | Fun_Leaf", "M_Rich | Fun_Root", "M_Rich | Fun_Rhizo",
                "M_Comp | Bac_Leaf", "M_Comp | Bac_Root", "M_Comp | Bac_Rhizo",
                "M_Comp | Fun_Leaf", "M_Comp | Fun_Root", "M_Comp | Fun_Rhizo")

row_levels_nl <- gsub(" \\| ", "\n", row_levels)

facet_order <- c("Plant performance", "Microbial richness", "Microbial community composition")

df <- df |>
  dplyr::mutate(
    Row_short  = forcats::fct_relevel(Row_short, row_levels),
    Row_short2 = factor(gsub(" \\| ", "\n", Row_short), levels = row_levels_nl),
    Group_full = forcats::fct_relevel(Group_full, facet_order),
    Factor_ord = factor(Factor, levels = rev(factor_levels)),
    value_lab   = sprintf("%.2f", Variance_explained_percent),
    value_lab_star = ifelse(Sig == "Yes", paste0(value_lab, "*"), value_lab),
    is_yes  = Sig == "Yes",
    is_marg = Sig == "Marginal")

df_yes   <- df |> dplyr::filter(is_yes)
df_marg  <- df |> dplyr::filter(is_marg)
df_other <- df |> dplyr::filter(!is_yes & !is_marg)

# -------------------------------------------------------------------
# 4) Plot heatmap of variance explained
# -------------------------------------------------------------------
p_var <- ggplot2::ggplot(df,
                         ggplot2::aes(x = Row_short2, y = Factor_ord, fill = Variance_explained_percent)) +
  ggplot2::geom_tile(color   = "grey40", linetype = "dotted", linewidth = 0.2) +
  ggplot2::geom_text(data = df_yes, ggplot2::aes(label = value_lab_star),
                     size = 3.6, fontface = "bold") +
  ggplot2::geom_text(data = df_marg, ggplot2::aes(label = value_lab),
                     size = 3.6, fontface = "italic") +
  ggplot2::geom_text(data = df_other, ggplot2::aes(label = value_lab), size = 3.6) +
  ggplot2::facet_grid(. ~ Group_full, scales = "free_x",
                      space  = "free_x", switch = "x") +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::scale_y_discrete(labels = setNames(factor_labels, factor_levels)) +
  ggplot2::scale_fill_gradient(name = "Variance explained (%)",
                               low  = "white", high = "royalblue") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x.top = ggplot2::element_text(
    face = "bold", size = 11, vjust = 0.3, color = "black"),
    axis.text.y = ggplot2::element_text(
      face = "bold", size = 13, color = "black"),
    axis.ticks.y = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.position  = "bottom",
    strip.text.x.bottom = ggplot2::element_text(
      face = "bold", size = 13, margin = ggplot2::margin(t = 4, b = 2)),
    panel.spacing = grid::unit(0.3, "lines"),
    panel.grid    = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title    = ggplot2::element_text(face = "bold", size = 12),
    legend.text     = ggplot2::element_text(size = 12))

p_var

# -------------------------------------------------------------------
# 5) Save the figure
# -------------------------------------------------------------------
ggplot2::ggsave(file.path(results_dir, "Variance_explained_heatmap.png"), ## Figure 5
                p_var, width = 10, height = 6, units = "in", dpi = 600, bg = "white")

