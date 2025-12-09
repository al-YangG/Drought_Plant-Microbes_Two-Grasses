##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/6b_taxonomic_top5_models.R
##
## NOTE:
##   This script fits models for the most abundant bacterial and fungal
##   taxa (top 5 phyla and orders) using per-sample relative abundances.
##   It uses beta GLMMs (glmmTMB) with:
##      Abundance ~ Compartment * Moisture_sc * Climate * Species
##      + (1 | Mesocosm)
##
## Tasks:
##   1) Load ps-Bac and ps-Fun phyloseq objects
##   2) Add MoistureLevel factor (M1–M6)
##   3) Build long data frames of relative abundance at Phylum and Order
##   4) Identify global top 5 phyla and orders (Bacteria/Fungi)
##   5) Fit beta GLMMs per taxon (Bacteria/Fungi; Phylum/Order)
##   6) Export model Anova tables to CSV

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root   <- here::here()
ps_dir <- file.path(root, "Data", "2_Phyloseq objects", "Output", "ps_new")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Load phyloseq objects
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(ps_dir, "ps-Bac.rds"))  # Bacteria
ps_Fun <- readRDS(file.path(ps_dir, "ps-Fun.rds"))  # Fungi

# -------------------------------------------------------------------
# 2) Add MoistureLevel (M1–M6)
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
# 3) Build per-sample relative abundance tables at a given rank
# -------------------------------------------------------------------
make_taxon_abund_df <- function(ps, rank) {
  ps_g <- phyloseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
  ps_rel <- phyloseq::transform_sample_counts(
    ps_g,
    function(x) if (sum(x) > 0) x / sum(x) else x
  )
  
  df <- phyloseq::psmelt(ps_rel) |>
    tibble::as_tibble()
  
  # Ensure standard types for predictors
  df |>
    dplyr::mutate(
      SampleID   = as.character(Sample),
      Mesocosm   = factor(Mesocosm),
      Species    = factor(Species),
      Climate    = factor(Climate),
      Compartment = factor(Compartment),
      Moisture   = as.numeric(Moisture),
      !!rank     := as.factor(.data[[rank]])
    )
}

# Bacteria
bac_phylum_df <- make_taxon_abund_df(ps_Bac, "phylum")
bac_order_df  <- make_taxon_abund_df(ps_Bac, "order")

# Fungi
fun_phylum_df <- make_taxon_abund_df(ps_Fun, "phylum")
fun_order_df  <- make_taxon_abund_df(ps_Fun, "order")

# -------------------------------------------------------------------
# 4) Identify global top 5 taxa (exclude unclassified)
# -------------------------------------------------------------------
top5_from_df <- function(df, tax_col,
                         exclude_regex = "incertae|unclassified|unknown|NA") {
  df |>
    dplyr::filter(
      !is.na(.data[[tax_col]]),
      !stringr::str_detect(
        .data[[tax_col]],
        stringr::regex(exclude_regex, ignore_case = TRUE)
      )
    ) |>
    dplyr::group_by(.data[[tax_col]]) |>
    dplyr::summarise(
      mean_abund = mean(Abundance, na.rm = TRUE),
      tot_abund  = sum(Abundance,  na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(mean_abund),
      dplyr::desc(tot_abund),
      .data[[tax_col]]
    ) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(.data[[tax_col]])
}

top5_bac_phylum <- top5_from_df(bac_phylum_df, "phylum")
top5_bac_order  <- top5_from_df(bac_order_df,  "order")

top5_fun_phylum <- top5_from_df(fun_phylum_df, "phylum")
top5_fun_order  <- top5_from_df(fun_order_df,  "order")

# Filter to top 5 taxa
bac_phylum_top5 <- bac_phylum_df |>
  dplyr::filter(phylum %in% top5_bac_phylum) |>
  droplevels()

bac_order_top5 <- bac_order_df |>
  dplyr::filter(order %in% top5_bac_order) |>
  droplevels()

fun_phylum_top5 <- fun_phylum_df |>
  dplyr::filter(phylum %in% top5_fun_phylum) |>
  droplevels()

fun_order_top5 <- fun_order_df |>
  dplyr::filter(order %in% top5_fun_order) |>
  droplevels()

# -------------------------------------------------------------------
# 5) Helper: beta GLMM per taxon
# -------------------------------------------------------------------
run_beta_models <- function(df, tax_col) {
  
  if (!"Mesocosm" %in% names(df)) {
    stop("Mesocosm column is required in the input data frame.")
  }
  
  df <- df |>
    dplyr::mutate(
      Moisture_sc = as.numeric(scale(Moisture)),
      Taxon       = .data[[tax_col]]
    )
  
  split_list <- split(df, df[[tax_col]])
  
  res_list <- lapply(split_list, function(sub) {
    
    # skip degenerate taxa (all zero/NA)
    if (all(is.na(sub$Abundance)) || max(sub$Abundance, na.rm = TRUE) <= 0) {
      return(NULL)
    }
    
    n <- nrow(sub)
    
    # Smithson & Verkuilen (2006) transformation to keep y in (0,1)
    sub$Abund01 <- (sub$Abundance * (n - 1) + 0.5) / n
    sub$Abund01[sub$Abund01 <= 0] <- 1e-6
    sub$Abund01[sub$Abund01 >= 1] <- 1 - 1e-6
    
    m <- glmmTMB::glmmTMB(
      Abund01 ~ Compartment * Moisture_sc * Climate * Species + (1 | Mesocosm),
      data   = sub,
      family = glmmTMB::beta_family(link = "logit")
    )
    
    a <- car::Anova(m, type = 3)
    
    out <- data.frame(
      Taxon   = unique(sub[[tax_col]]),
      Effect  = rownames(a),
      Chisq   = a$Chisq,
      Df      = a$Df,
      p_value = a$`Pr(>Chisq)`,
      row.names = NULL
    )
    
    out
  })
  
  dplyr::bind_rows(res_list)
}

# -------------------------------------------------------------------
# 6) Run models and export results
# -------------------------------------------------------------------
## Bacteria – Phylum
res_bac_phylum_beta <- run_beta_models(bac_phylum_top5, "phylum")
write.csv(res_bac_phylum_beta,
          file.path(results_dir, "Top5_Bacteria_Phylum_beta.csv"), ## Table S10a
          row.names = FALSE)

## Bacteria – Order
res_bac_order_beta <- run_beta_models(bac_order_top5, "order")
write.csv(res_bac_order_beta,
          file.path(results_dir, "Top5_Bacteria_Order_beta.csv"), ## Table S10c
          row.names = FALSE)

## Fungi – Phylum
res_fun_phylum_beta <- run_beta_models(fun_phylum_top5, "phylum")
write.csv(res_fun_phylum_beta,
          file.path(results_dir, "Top5_Fungi_Phylum_beta.csv"), ## Table S10b
          row.names = FALSE)

## Fungi – Order
res_fun_order_beta <- run_beta_models(fun_order_top5, "order")
write.csv(res_fun_order_beta,
          file.path(results_dir, "Top5_Fungi_Order_beta.csv"), ## Table S10d
          row.names = FALSE)
