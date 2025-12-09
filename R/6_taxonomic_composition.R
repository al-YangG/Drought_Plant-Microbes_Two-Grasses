##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/6_taxonomic_composition.R
##
## NOTE:
##   This script describes taxonomic composition (relative abundances)
##   of bacteria and fungi, and identifies the top 5 phyla and orders
##   across compartments, moisture levels, species and climate.
##
## Tasks:
##   1) Load ps-Bac and ps-Fun phyloseq objects
##   2) Add MoistureLevel factor (M1–M6) to sample_data
##   3) Helper: rel_topN_by_group() for top-N relative abundances
##   4) Compute global top-5 coverage across taxonomic ranks
##   5) Compute top-5 phyla/orders by group (Compartment, MoistureLevel,
##      Species, Climate) for Bacteria and Fungi
##   6) Save grouped top-5 tables
##   7) Plot stacked-bar compositions for top 5 phyla and orders
##      (Bacteria + Fungi) and save combined figure

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
ps_dir      <- file.path(root, "Data", "2_Phyloseq objects", "Output", "ps_new")
data_dir    <- file.path(root, "Data", "6_Taxonomic composition")
results_dir <- file.path(root, "results")

dir.create(data_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Load phyloseq objects (domain-level)
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(ps_dir, "ps-Bac.rds"))  # Bacteria
ps_Fun <- readRDS(file.path(ps_dir, "ps-Fun.rds"))  # Fungi

# -------------------------------------------------------------------
# 2) Add MoistureLevel (M1–M6) to sample_data
# -------------------------------------------------------------------
add_M_levels <- function(ps) {
  sd <- as.data.frame(phyloseq::sample_data(ps))
  sd$Moisture <- as.numeric(sd$Moisture)
  
  sd$MoistureLevel <- dplyr::case_when(
    sd$Moisture %in% 0              ~ "M1",
    sd$Moisture %in% c(300, 400)    ~ "M2",
    sd$Moisture %in% c(600, 700)    ~ "M3",
    sd$Moisture %in% c(900, 1000)   ~ "M4",
    sd$Moisture %in% c(1200, 1300)  ~ "M5",
    sd$Moisture %in% c(1500, 1600)  ~ "M6",
    TRUE                            ~ NA_character_
  )
  
  sd$MoistureLevel <- factor(sd$MoistureLevel, levels = paste0("M", 1:6))
  phyloseq::sample_data(ps)$MoistureLevel <- sd$MoistureLevel
  ps
}

ps_Bac <- add_M_levels(ps_Bac)
ps_Fun <- add_M_levels(ps_Fun)

# -------------------------------------------------------------------
# 3) Helper: Top-N relative abundances by group
# -------------------------------------------------------------------
rel_topN_by_group <- function(ps,
                              rank          = "phylum",
                              group_vars    = c("Compartment"),
                              top_n         = 5,
                              other_label   = "Other",
                              exclude_regex = "incertae|unclassified|unknown|NA") {
  
  stopifnot(inherits(ps, "phyloseq"))
  
  # Aggregate to chosen rank and convert to relative abundance
  ps_g   <- phyloseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
  ps_rel <- phyloseq::transform_sample_counts(ps_g,
                                              function(x) if (sum(x) > 0) x / sum(x) else x)
  
  df <- phyloseq::psmelt(ps_rel) |>
    tibble::as_tibble()
  
  # keep a character version of the rank column
  df <- df |>
    dplyr::mutate(.taxon_chr = as.character(.data[[rank]]))
  
  # 1) Determine the global top-N taxa (exclude unclassified from the choice)
  df_for_top <- df |>
    dplyr::filter(!is.na(.taxon_chr),
                  !stringr::str_detect(
                    .taxon_chr,
                    stringr::regex(exclude_regex, ignore_case = TRUE)))
  
  top_set <- df_for_top |>
    dplyr::group_by(.taxon_chr) |>
    dplyr::summarise(global_mean = mean(Abundance, na.rm = TRUE),
                     global_sum  = sum(Abundance,  na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(global_mean), dplyr::desc(global_sum),
                   .taxon_chr) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(.taxon_chr)
  
  # 2) Lump all non-top taxa to "Other"
  df_lumped <- df |>
    dplyr::mutate(Taxon = dplyr::if_else(.taxon_chr %in% top_set, .taxon_chr, other_label))
  
  # 3) Summarise mean relative abundance and SE
  if (is.null(group_vars) || length(group_vars) == 0) {
    
    out <- df_lumped |>
      dplyr::group_by(Taxon) |>
      dplyr::summarise(mean_RA = mean(Abundance, na.rm = TRUE),
                       se_RA   = stats::sd(Abundance, na.rm = TRUE) / sqrt(dplyr::n()),
                       .groups = "drop")
    
    # Re-normalize so totals sum exactly to 1
    denom <- sum(out$mean_RA, na.rm = TRUE)
    if (is.finite(denom) && denom > 0) {
      out <- out |>
        dplyr::mutate(mean_RA = mean_RA / denom, se_RA   = se_RA   / denom)
    }
    
  } else {
    
    out <- df_lumped |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)), Taxon) |>
      dplyr::summarise(
        mean_RA = mean(Abundance, na.rm = TRUE),
        se_RA   = stats::sd(Abundance, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups = "drop_last"
      ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::mutate(
        .denom  = sum(mean_RA, na.rm = TRUE),
        mean_RA = dplyr::if_else(.denom > 0, mean_RA / .denom, mean_RA),
        se_RA   = dplyr::if_else(.denom > 0, se_RA   / .denom, se_RA)
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.denom)
  }
  
  out
}

# -------------------------------------------------------------------
# 4) Global top-5 coverage by rank (Bacteria & Fungi)
# -------------------------------------------------------------------
# Ranks to evaluate
ranks_to_check <- c("phylum", "class", "order", "family", "genus")

# Bacteria
bac_top5_global <- lapply(ranks_to_check, function(rk) {
  tmp <- rel_topN_by_group(ps_Bac, rank = rk, group_vars = NULL, top_n = 5)
  tmp$RankLevel <- tools::toTitleCase(rk)
  tmp
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Domain = "Bacteria", .before = 1)

bac_coverage <- bac_top5_global |>
  dplyr::group_by(Domain, RankLevel) |>
  dplyr::summarise(Top5_Coverage = sum(mean_RA, na.rm = TRUE), .groups = "drop")

# Fungi
fun_top5_global <- lapply(ranks_to_check, function(rk) {
  tmp <- rel_topN_by_group(ps_Fun, rank = rk, group_vars = NULL, top_n = 5)
  tmp$RankLevel <- tools::toTitleCase(rk)
  tmp
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Domain = "Fungi", .before = 1)

fun_coverage <- fun_top5_global |>
  dplyr::group_by(Domain, RankLevel) |>
  dplyr::summarise(Top5_Coverage = sum(mean_RA, na.rm = TRUE), .groups = "drop")

coverage_all <- dplyr::bind_rows(bac_coverage, fun_coverage)

write.csv(coverage_all, file.path(results_dir, "Taxonomic_Top5_Coverage.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 5) Grouped top-5 tables (Phylum + Order; Bac + Fun)
# -------------------------------------------------------------------
group_vars <- c("MoistureLevel", "Compartment", "Species", "Climate")

## Bacteria
bac_phylum_top5 <- rel_topN_by_group(ps_Bac,
                                     rank       = "phylum",
                                     group_vars = group_vars,
                                     top_n      = 5) |>
  dplyr::mutate(Domain = "Bacteria",
                RankLevel = "Phylum",
                .before = 1)

bac_order_top5 <- rel_topN_by_group(ps_Bac,
                                    rank       = "order",
                                    group_vars = group_vars,
                                    top_n      = 5) |>
  dplyr::mutate(Domain = "Bacteria", RankLevel = "Order", .before = 1)

bac_top5_group <- dplyr::bind_rows(bac_phylum_top5, bac_order_top5)

write.csv(bac_top5_group, file.path(results_dir, "Taxonomic_Bacteria_Top5_byGroup.csv"),
          row.names = FALSE)

## Fungi
fun_phylum_top5 <- rel_topN_by_group(ps_Fun,
                                     rank       = "phylum",
                                     group_vars = group_vars,
                                     top_n      = 5) |>
  dplyr::mutate(Domain = "Fungi", RankLevel = "Phylum", .before = 1)

fun_order_top5 <- rel_topN_by_group(ps_Fun,
                                    rank       = "order",
                                    group_vars = group_vars,
                                    top_n      = 5) |>
  dplyr::mutate(Domain = "Fungi", RankLevel = "Order", .before = 1)

fun_top5_group <- dplyr::bind_rows(fun_phylum_top5, fun_order_top5)

write.csv(fun_top5_group, file.path(results_dir, "Taxonomic_Fungi_Top5_byGroup.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 6) Palettes and themes for plots
# -------------------------------------------------------------------
pal_bac_phylum <- c(
  "#1B9E77",  # teal
  "#E7298A",  # pink-magenta
  "#7570B3",  # purple
  "#D95F02",  # orange
  "#66A61E",  # green
  "#B3B3B3"   # other
)

pal_bac_order <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # violet
  "#FF7F00",  # orange
  "#B3B3B3"   # other
)

pal_fun_phylum <- c(
  "#A6761D",  # brown-gold
  "#E6AB02",  # mustard
  "#7570B3",  # violet
  "#66A61E",  # green
  "#E7298A",  # pink
  "#B3B3B3"   # other
)

pal_fun_order <- c(
  "#D73027",  # crimson
  "#FC8D59",  # salmon
  "#4575B4",  # blue
  "#91BFDB",  # light blue
  "#F46D43",  # coral-orange
  "#B3B3B3"   # other
)

# A compact facet theme for horizontal barplots
tax_theme <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.line         = ggplot2::element_blank(),
    axis.ticks.y      = ggplot2::element_blank(),
    axis.text.x       = ggplot2::element_text(face = "bold", size = 10, colour = "black"),
    axis.text.y       = ggplot2::element_text(face = "bold", size = 10, colour = "black"),
    axis.title        = ggplot2::element_text(face = "bold", size = 12, colour = "black"),
    legend.position   = "bottom",
    legend.text       = ggplot2::element_text(face = "bold", size = 10),
    legend.title      = ggplot2::element_text(face = "bold", size = 10),
    strip.text        = ggplot2::element_text(face = "bold", size = 11),
    strip.background  = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.6),
    panel.background  = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.6),
    panel.spacing     = grid::unit(0.5, "lines"),
    plot.title        = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
  )

# -------------------------------------------------------------------
# 7) Helper to build horizontal stacked barplots for top-5
# -------------------------------------------------------------------
build_top5_barplot <- function(df_top5, rank_level, palette, rank_label) {
  
  df_rank <- df_top5 |>
    dplyr::filter(RankLevel == rank_level)
  
  # global taxon order (Other last)
  tax_order <- df_rank |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(m = mean(mean_RA, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(m)) |>
    dplyr::pull(Taxon)
  
  tax_order <- c(setdiff(tax_order, "Other"), "Other")
  
  # 1) Compartment
  comp <- df_rank |>
    dplyr::group_by(Compartment, Taxon) |>
    dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      Compartment = dplyr::recode(
        Compartment,
        "Leaf" = "Leaf endosphere",
        "Root" = "Root endosphere",
        "Soil" = "Rhizosphere"
      ),
      Compartment = factor(
        Compartment,
        levels = c("Rhizosphere", "Root endosphere", "Leaf endosphere")
      )
    ) |>
    dplyr::group_by(Compartment) |>
    dplyr::mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(Taxon = factor(Taxon, levels = tax_order))
  
  # 2) Moisture level
  moist <- df_rank |>
    dplyr::group_by(MoistureLevel, Taxon) |>
    dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      MoistureLevel = factor(MoistureLevel, levels = paste0("M", 6:1))
    ) |>
    dplyr::group_by(MoistureLevel) |>
    dplyr::mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(Taxon = factor(Taxon, levels = tax_order))
  
  # 3) Species
  sp <- df_rank |>
    dplyr::group_by(Species, Taxon) |>
    dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      Species = factor(Species, levels = c("Lolium", "Festuca"))
    ) |>
    dplyr::group_by(Species) |>
    dplyr::mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(Taxon = factor(Taxon, levels = tax_order))
  
  # 4) Climate
  clim <- df_rank |>
    dplyr::group_by(Climate, Taxon) |>
    dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
      Climate = factor(Climate, levels = c("Current", "Future"))
    ) |>
    dplyr::group_by(Climate) |>
    dplyr::mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(Taxon = factor(Taxon, levels = tax_order))
  
  # combine with group type and "X" variable
  comb <- dplyr::bind_rows(
    comp  |> dplyr::transmute(GroupType = "Compartment",    X = Compartment,   Taxon, mean_RA),
    moist |> dplyr::transmute(GroupType = "Moisture level", X = MoistureLevel, Taxon, mean_RA),
    sp    |> dplyr::transmute(GroupType = "Species",        X = Species,       Taxon, mean_RA),
    clim  |> dplyr::transmute(GroupType = "Climate",        X = Climate,       Taxon, mean_RA)
  ) |>
    dplyr::mutate(
      GroupType = factor(
        GroupType,
        levels = c("Compartment", "Moisture level", "Species", "Climate")
      ),
      Taxon = factor(Taxon, levels = rev(tax_order))
    )
  
  # ensure X is ordered inside each GroupType (reverse for nicer vertical layout)
  comb <- comb |>
    dplyr::group_by(GroupType) |>
    dplyr::mutate(
      X = factor(X, levels = rev(unique(as.character(X))))
    ) |>
    dplyr::ungroup()
  
  ggplot2::ggplot(comb, ggplot2::aes(x = mean_RA, y = X, fill = Taxon)) +
    ggplot2::geom_col(
      width    = 0.7,
      position = "stack",
      colour   = NA
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE, guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x    = "Mean relative abundance (%)",
      y    = NULL,
      fill = rank_label
    ) +
    ggplot2::facet_grid(
      rows   = ggplot2::vars(GroupType),
      scales = "free_y",
      space  = "free_y"
    ) +
    tax_theme
}

# -------------------------------------------------------------------
# 8) Build plots and save combined figure
# -------------------------------------------------------------------
# Re-use grouped tables already in memory
bac_top5_group <- bac_top5_group
fun_top5_group <- fun_top5_group

p_bac_phylum <- build_top5_barplot(df_top5    = bac_top5_group,
                                   rank_level = "Phylum",
                                   palette    = pal_bac_phylum,
                                   rank_label = "Bacterial phylum")

p_fun_phylum <- build_top5_barplot(df_top5    = fun_top5_group,
                                   rank_level = "Phylum",
                                   palette    = pal_fun_phylum,
                                   rank_label = "Fungal phylum")

p_bac_order <- build_top5_barplot(df_top5    = bac_top5_group,
                                  rank_level = "Order",
                                  palette    = pal_bac_order,
                                  rank_label = "Bacterial order")

p_fun_order <- build_top5_barplot(df_top5    = fun_top5_group,
                                  rank_level = "Order",
                                  palette    = pal_fun_order,
                                  rank_label = "Fungal order")

# combined 2 × 2 panel
p_tax_all <- ggpubr::ggarrange(p_bac_phylum, p_fun_phylum,
                               p_bac_order,  p_fun_order,
                               nrow   = 2, ncol   = 2,
                               labels = c("(a)", "(b)", "(c)", "(d)"))

ggplot2::ggsave(file.path(results_dir, "Taxonomic_Top5_Phyla_Orders.png"), ## Figure S11
                p_tax_all, width = 13, height = 15, units = "in",
                bg = "white", dpi = 600)

