##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/2b_qc_venn_rarefaction.R
##
## Tasks:
##  1) Venn diagrams (Bacteria / Fungi by Compartment, Species, Climate)
##  2) Sequencing depth QC:
##       - depth summary
##       - survival vs thresholds
##       - rarefaction curves (Figure S2)

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()   # ← loads phyloseq, dplyr, tidyr, MicEco, vegan, ggplot2, etc.

root     <- here::here()
data_dir <- file.path(root, "Data", "2_Phyloseq objects")
ps_dir   <- file.path(data_dir, "Output", "ps_new")
res_dir  <- file.path(root, "results")

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# color palettes (from utils_theme.R)
col_species <- species_cols
col_climate <- climate_cols
col_comp    <- c("#afffa8", "#ffd966", "#f2f2f2")

# -------------------------------------------------------------------
# 1) Load phyloseq objects
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(ps_dir, "ps-Bac.rds"))
ps_Fun <- readRDS(file.path(ps_dir, "ps-Fun.rds"))

ps_BL <- readRDS(file.path(ps_dir, "ps-BL.rds"))
ps_BR <- readRDS(file.path(ps_dir, "ps-BR.rds"))
ps_BS <- readRDS(file.path(ps_dir, "ps-BS.rds"))
ps_FL <- readRDS(file.path(ps_dir, "ps-FL.rds"))
ps_FR <- readRDS(file.path(ps_dir, "ps-FR.rds"))
ps_FS <- readRDS(file.path(ps_dir, "ps-FS.rds"))

ps_list <- list(
  "Bacteria - Leaf endosphere"  = ps_BL,
  "Bacteria - Root endosphere"  = ps_BR,
  "Bacteria - Rhizosphere"      = ps_BS,
  "Fungi - Leaf endosphere"     = ps_FL,
  "Fungi - Root endosphere"     = ps_FR,
  "Fungi - Rhizosphere"         = ps_FS
)

# -------------------------------------------------------------------
# 2) Venn diagrams
# -------------------------------------------------------------------

## Compartment
vb1 <- ps_venn(ps_Bac, group = "Compartment", fill = col_comp)
vf1 <- ps_venn(ps_Fun, group = "Compartment", fill = col_comp)

## Species
vb2 <- ps_venn(ps_Bac, group = "Species",  fill = col_species)
vf2 <- ps_venn(ps_Fun, group = "Species",  fill = col_species)

## Climate
vb3 <- ps_venn(ps_Bac, group = "Climate",  fill = col_climate)
vf3 <- ps_venn(ps_Fun, group = "Climate",  fill = col_climate)

# arrange panels
v_all <- ggpubr::ggarrange(
  ggpubr::ggarrange(vb1, vb3, vb2, ncol = 1, labels = c("(a)", "(c)", "(e)")),
  ggpubr::ggarrange(vf1, vf3, vf2, ncol = 1, labels = c("(b)", "(d)", "(f)")),
  nrow = 1, widths = c(1, 1)
)

ggsave(file.path(res_dir, "Venn_Bacteria_Fungi.png"), ## Figure S4
       v_all, bg = "white", width = 8, height = 10, dpi = 600)

# -------------------------------------------------------------------
# 3) Rarefaction QC
# -------------------------------------------------------------------

## Helpers
otu_to_matrix <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

depths_df_all <- function(ps_list) {
  bind_rows(lapply(names(ps_list), function(nm) {
    tibble(Dataset  = nm, SampleID = sample_names(ps_list[[nm]]),
           Depth    = sample_sums(ps_list[[nm]]))
  }))
}

# Depth summary
depths_all <- depths_df_all(ps_list)

write.csv(depths_all %>%
            group_by(Dataset) %>% 
            summarise(min = min(Depth), q1 = quantile(Depth, .25),
                      median = median(Depth), mean = mean(Depth),
                      q3 = quantile(Depth, .75), max = max(Depth)),
          file.path(res_dir, "Depth_summary_all_compartments.csv"), ## Table S2
          row.names = FALSE)

# Survival vs proposed thresholds
thresholds <- c(1000, 2500, 5000, 10000, 15000)
survival_tab <- depths_all %>%
  tidyr::crossing(Threshold = thresholds) %>%
  mutate(Keep = Depth >= Threshold) %>%
  group_by(Dataset, Threshold) %>%
  summarise(n_total = n(), n_kept = sum(Keep),
            kept_pct = round(100 * n_kept / n_total, 1), .groups = "drop")

write.csv(survival_tab,
          file.path(res_dir, "Survival_by_threshold_all_compartments.csv"), ## Table S3
          row.names = FALSE)

# Rarefaction curves
rare_curve_one_sample <- function(counts, sample_id, n_points = 40) {
  lib <- sum(counts)
  if (lib < 10) return(tibble())
  depths <- unique(round(exp(seq(log(10), log(lib), length.out = n_points))))
  S <- vegan::rarefy(matrix(counts, nrow = 1), sample = depths)
  tibble(SampleID = sample_id, Depth = depths, ObservedExp = as.numeric(S))
}

compute_rare_df_ps <- function(ps, dataset_name) {
  mat <- otu_to_matrix(ps)
  ids <- rownames(mat)
  purrr::map2_dfr(asplit(mat, 1), ids,
                  ~ rare_curve_one_sample(as.numeric(.x), .y)) %>%
    mutate(Dataset = dataset_name)
}

rare_df <- purrr::imap_dfr(ps_list, compute_rare_df_ps)

# median rarefaction curve per dataset
median_df <- rare_df %>%
  group_by(Dataset) %>%
  group_modify(~ {
    upper <- max(10, quantile(.x$Depth, .95, na.rm = TRUE))
    grid <- unique(round(exp(seq(log(10), log(upper), length.out = 60))))
    by_samp <- group_split(.x, .x$SampleID)
    interp <- lapply(by_samp, function(d) {
      d <- arrange(d, Depth)
      if (nrow(d) >= 2) {
        out <- approx(x = d$Depth, y = d$ObservedExp,
                      xout = grid, method = "linear", rule = 2)
        tibble(Depth = out$x, ObservedExp = out$y)
      } else {
        tibble(Depth = grid,
               ObservedExp = ifelse(grid <= max(d$Depth), d$ObservedExp[1], NA_real_))
      }
    })
    bind_rows(interp) %>%
      group_by(Depth) %>%
      summarise(MedianObserved = median(ObservedExp, na.rm = TRUE), .groups = "drop")
  })

cutoffs <- depths_all %>%
  group_by(Dataset) %>%
  summarise(cutoff = min(Depth), .groups = "drop")

order_vec <- names(ps_list)
rare_df$Dataset   <- factor(rare_df$Dataset,   levels = order_vec)
median_df$Dataset <- factor(median_df$Dataset, levels = order_vec)
cutoffs$Dataset   <- factor(cutoffs$Dataset,   levels = order_vec)

p_rare <- ggplot(rare_df, aes(Depth, ObservedExp, group = SampleID)) +
  geom_path(alpha = 0.5, linewidth = 0.3, color = "grey40") +
  geom_path(data = median_df,
            aes(Depth, MedianObserved),
            color = "red", linewidth = 1.1) +
  geom_vline(data = cutoffs,
             aes(xintercept = cutoff),
             linetype = "dashed", color = "blue") +
  facet_wrap(~ Dataset, scales = "free_y", ncol = 3) +
  labs(
    title = "Rarefaction curves by compartment",
    x = "Sequencing depth (reads)",
    y = "Observed ASVs (expected)"
  ) +
  project_theme_classic()

ggsave(file.path(res_dir, "Rarefaction_curves.png"), ## Figure S3
       p_rare, bg = "white", width = 10, height = 6, dpi = 600)

