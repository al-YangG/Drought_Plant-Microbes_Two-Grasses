##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/2_phyloseq_objects.R
##
## Tasks:
##  1) Construct phyloseq objects (BL, BR, BS, FL, FR, FS) from CSVs
##  2) Merge to domain-level ps objects (ps-Bac, ps-Fun)
##  3) Rarefy BR, BS, FL, FR, FS using utils_phyloseq::rarefy_ps_safe()
##  4) Compute alpha-diversity metrics with utils_phyloseq::compute_alpha_div()
##  5) Export per-compartment alpha tables + combined alpha table

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")
source("R/utils_phyloseq.R")

load_project_packages()

root      <- here::here()
data_dir  <- file.path(root, "Data", "2_Phyloseq objects")
input_dir <- file.path(data_dir, "Input")
asv_dir   <- file.path(input_dir, "asv")
tax_dir   <- file.path(input_dir, "taxa")
env_dir   <- file.path(input_dir, "env")

ps_dir    <- file.path(data_dir, "Output", "ps_new")
alpha_dir <- file.path(root, "results")

dir.create(ps_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(alpha_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Construct phyloseq objects for each compartment
#    (uses build_phyloseq_from_csv() from utils_phyloseq.R)
# -------------------------------------------------------------------

## Bacteria
ps_bl <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-BL.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-BL.csv"),
  meta_fp = file.path(env_dir,  "env-BL.csv")
)

ps_br <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-BR.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-BR.csv"),
  meta_fp = file.path(env_dir,  "env-BR.csv")
)

ps_bs <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-BS.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-BS.csv"),
  meta_fp = file.path(env_dir,  "env-BS.csv")
)

## Fungi
ps_fl <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-FL.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-FL.csv"),
  meta_fp = file.path(env_dir,  "env-FL.csv")
)

ps_fr <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-FR.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-FR.csv"),
  meta_fp = file.path(env_dir,  "env-FR.csv")
)

ps_fs <- build_phyloseq_from_csv(
  asv_fp  = file.path(asv_dir,  "asv-FS.csv"),
  taxa_fp = file.path(tax_dir, "taxonomy-FS.csv"),
  meta_fp = file.path(env_dir,  "env-FS.csv")
)

## Save individual ps objects (compartment-level)
saveRDS(ps_bl, file.path(ps_dir, "ps-BL.rds"))
saveRDS(ps_br, file.path(ps_dir, "ps-BR.rds"))
saveRDS(ps_bs, file.path(ps_dir, "ps-BS.rds"))
saveRDS(ps_fl, file.path(ps_dir, "ps-FL.rds"))
saveRDS(ps_fr, file.path(ps_dir, "ps-FR.rds"))
saveRDS(ps_fs, file.path(ps_dir, "ps-FS.rds"))

# -------------------------------------------------------------------
# 2) Merge to domain-level ps objects (ps-Bac, ps-Fun)
# -------------------------------------------------------------------
ps_Bac <- phyloseq::merge_phyloseq(ps_bl, ps_br, ps_bs)
ps_Fun <- phyloseq::merge_phyloseq(ps_fl, ps_fr, ps_fs)

saveRDS(ps_Bac, file.path(ps_dir, "ps-Bac.rds"))
saveRDS(ps_Fun, file.path(ps_dir, "ps-Fun.rds"))

# -------------------------------------------------------------------
# 3) Rarefaction per compartment (using rarefy_ps_safe)
#     We follow your original choice: no alpha diversity for BL.
# -------------------------------------------------------------------
# Reload from disk to keep I/O boundary clear
ps_br <- readRDS(file.path(ps_dir, "ps-BR.rds"))
ps_bs <- readRDS(file.path(ps_dir, "ps-BS.rds"))
ps_fl <- readRDS(file.path(ps_dir, "ps-FL.rds"))
ps_fr <- readRDS(file.path(ps_dir, "ps-FR.rds"))
ps_fs <- readRDS(file.path(ps_dir, "ps-FS.rds"))

ps_BR_r <- rarefy_ps_safe(ps_br, seed = 100)
ps_BS_r <- rarefy_ps_safe(ps_bs, seed = 100)
ps_FL_r <- rarefy_ps_safe(ps_fl, seed = 100)
ps_FR_r <- rarefy_ps_safe(ps_fr, seed = 100)
ps_FS_r <- rarefy_ps_safe(ps_fs, seed = 100)

# (optional) save rarefied ps objects
saveRDS(ps_BR_r, file.path(ps_dir, "ps-BR_rarefied.rds"))
saveRDS(ps_BS_r, file.path(ps_dir, "ps-BS_rarefied.rds"))
saveRDS(ps_FL_r, file.path(ps_dir, "ps-FL_rarefied.rds"))
saveRDS(ps_FR_r, file.path(ps_dir, "ps-FR_rarefied.rds"))
saveRDS(ps_FS_r, file.path(ps_dir, "ps-FS_rarefied.rds"))

# -------------------------------------------------------------------
# 4) Alpha-diversity metrics (compute_alpha_div + join metadata)
# -------------------------------------------------------------------
alpha_from_ps_with_meta <- function(ps) {
  # diversity matrix (rows = samples)
  div_mat <- compute_alpha_div(ps)           # from utils_phyloseq.R
  div_df  <- as.data.frame(div_mat)
  div_df$SampleID <- rownames(div_df)
  
  # sample metadata
  env_df <- as.data.frame(phyloseq::sample_data(ps))
  env_df$SampleID <- rownames(env_df)
  
  dplyr::left_join(env_df, div_df, by = "SampleID")
}

alpha_BR <- alpha_from_ps_with_meta(ps_BR_r)  # Bacteria - Root
alpha_BS <- alpha_from_ps_with_meta(ps_BS_r)  # Bacteria - Rhizosphere
alpha_FL <- alpha_from_ps_with_meta(ps_FL_r)  # Fungi - Leaf
alpha_FR <- alpha_from_ps_with_meta(ps_FR_r)  # Fungi - Root
alpha_FS <- alpha_from_ps_with_meta(ps_FS_r)  # Fungi - Rhizosphere

# -------------------------------------------------------------------
# 5) Export per-compartment alpha tables
# -------------------------------------------------------------------
write.csv(alpha_BR,
          file.path(alpha_dir, "Alpha-BR.csv"),
          row.names = FALSE)

write.csv(alpha_BS,
          file.path(alpha_dir, "Alpha-BS.csv"),
          row.names = FALSE)

write.csv(alpha_FL,
          file.path(alpha_dir, "Alpha-FL.csv"),
          row.names = FALSE)

write.csv(alpha_FR,
          file.path(alpha_dir, "Alpha-FR.csv"),
          row.names = FALSE)

write.csv(alpha_FS,
          file.path(alpha_dir, "Alpha-FS.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 6) Combined alpha-diversity table across compartments
# -------------------------------------------------------------------
alpha_all <- dplyr::bind_rows(
  dplyr::mutate(alpha_BR, Dataset = "BR (Bacteria-Root)"),
  dplyr::mutate(alpha_BS, Dataset = "BS (Bacteria-Rhizosphere)"),
  dplyr::mutate(alpha_FL, Dataset = "FL (Fungi-Leaf)"),
  dplyr::mutate(alpha_FR, Dataset = "FR (Fungi-Root)"),
  dplyr::mutate(alpha_FS, Dataset = "FS (Fungi-Rhizosphere)")
) |>
  dplyr::relocate(Dataset, SampleID)

write.csv(alpha_all, file.path(alpha_dir, "Alpha-all_compartments.csv"), row.names = FALSE)


