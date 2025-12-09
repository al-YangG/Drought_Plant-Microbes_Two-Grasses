##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/8_functional_pathways_Tax4Fun2.R
##
## Tasks:
##  1) (Optional) Run Tax4Fun2 functional prediction for bacterial 16S data.
##  2) Merge KEGG pathway predictions with sample metadata and export
##     "Pathway with metadata_cleaned.csv".
##  3) From the cleaned table, extract top 15 pathways by:
##        - Plant species (Festuca, Lolium)
##        - Compartment (Leaf / Root / Rhizosphere)
##        - Moisture (6 drought classes M1–M6)
##  4) Plot stacked barplots and export combined figure.
##
## Note:
##  - Adjust paths to fasta / OTU table / reference DB / metadata as needed.
##  - If you already have "Pathway with metadata_cleaned.csv", you can
##    skip Section 2 and start from Section 3.
## -------------------------------------------------------------------

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
func_dir    <- file.path(root, "Data", "8_Functional group")
tax4fun_out <- file.path(func_dir, "Tax4Fun2_output_folder")
results_dir <- file.path(root, "results")

dir.create(tax4fun_out, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Tax4Fun2 functional prediction (optional, run once)
# -------------------------------------------------------------------
## Requires:
##  - Bac-seqs.fasta
##  - asv-filterted.tsv          (ASV table; samples in columns)
##  - Tax4Fun2_ReferenceData_v2  (downloaded reference DB)
##
## If predictions already done and "pathway_prediction.csv" exists,
## you can skip this section.

## Load Tax4Fun2 explicitly (not in utils_packages)
if (!requireNamespace("Tax4Fun2", quietly = TRUE)) {
  stop("Package 'Tax4Fun2' not installed. Install it or comment out this section.")
}
library(Tax4Fun2)

# Input files (adjust if needed)
query_otu_seq   <- file.path(func_dir, "Bac-seqs.fasta")
query_otu_table <- file.path(func_dir, "asv-filterted.tsv")
pwd_ref_data    <- file.path(func_dir, "Tax4Fun2_ReferenceData_v2")

# Settings
norm_by_cn  <- TRUE   # normalize_by_copy_number
norm_path   <- TRUE   # normalize_pathways
iden        <- 0.97   # min identity to reference
num_threads <- 6      # CPU cores

# 1.1 runRefBlast
runRefBlast(path_to_otus = query_otu_seq,
            path_to_reference_data = pwd_ref_data,
            path_to_temp_folder = tax4fun_out,
            database_mode = "Ref99NR",
            use_force = TRUE,
            num_threads = num_threads)

# 1.2 functional prediction
makeFunctionalPrediction(path_to_otu_table = query_otu_table,
                         path_to_reference_data = pwd_ref_data,
                         path_to_temp_folder = tax4fun_out,
                         database_mode = "Ref99NR",
                         normalize_by_copy_number = norm_by_cn,
                         min_identity_to_reference = iden,
                         normalize_pathways = norm_path)

# -------------------------------------------------------------------
# 2) Build “Pathway with metadata_cleaned.csv”
# -------------------------------------------------------------------
## If this file already exists (in Data/8_Functional group), you can
## skip to Section 3.

# Predicted pathway table from Tax4Fun2
pathway_fp <- file.path(tax4fun_out, "pathway_prediction.csv")
pathway    <- read.csv(pathway_fp, check.names = FALSE)

# Reshape to long format: one row per pathway × sample
pathway_long <- pathway |>
  tidyr::pivot_longer(cols = -c(pathway, level1, level2, level3),
                      names_to  = "SampleID", values_to = "Abundance")

# Metadata (adjust path if your metadata is stored elsewhere)
meta_fp <- file.path(root, "Data", "metadata", "Metadata_all_samples.csv")
metadata <- read.csv(meta_fp, check.names = FALSE)

# Merge with metadata
merged_df <- pathway_long |>
  dplyr::left_join(metadata, by = "SampleID")

# Remove pathways lacking level2/level3 annotation
cleaned_pathways <- merged_df |>
  dplyr::filter(!is.na(level2), level2 != "",
                !is.na(level3), level3 != "") |>
  dplyr::mutate(RelativeAbundance = Abundance * 100)

# Save cleaned table (keep original filename with spaces)
cleaned_fp <- file.path(func_dir, "Pathway with metadata_cleaned.csv")
write.csv(cleaned_pathways, cleaned_fp, row.names = FALSE)

# -------------------------------------------------------------------
# 3) Read cleaned pathway table
# -------------------------------------------------------------------
pathway_cleaned <- read.csv(cleaned_fp, check.names = FALSE)
str(pathway_cleaned)

# -------------------------------------------------------------------
# 4) Helper: top-N pathways per group
# -------------------------------------------------------------------
# generic helper: compute mean RelativeAbundance and take top N per group
top_n_pathways <- function(df, group_var, n = 15) {
  group_sym <- rlang::sym(group_var)
  
  df |>
    dplyr::group_by(!!group_sym, pathway, level1, level2, level3) |>
    dplyr::summarise(Average_Abundance = mean(RelativeAbundance, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::group_by(!!group_sym) |>
    dplyr::slice_max(order_by  = Average_Abundance, n = n, with_ties = FALSE) |>
    dplyr::ungroup()
}

# -------------------------------------------------------------------
# 5) Top 15 by plant species (Festuca vs Lolium)
# -------------------------------------------------------------------
top_species <- pathway_cleaned |>
  dplyr::filter(!is.na(Species)) |>
  top_n_pathways(group_var = "Species", n = 15) |>
  dplyr::mutate(Species = factor(Species, levels = c("Festuca", "Lolium")))

# rename for plotting
top_s <- top_species

# -------------------------------------------------------------------
# 6) Top 15 by compartment (Leaf / Root / Rhizosphere)
# -------------------------------------------------------------------
# First harmonise compartment names if needed
top_compartment <- pathway_cleaned |>
  dplyr::mutate(Compartment = dplyr::case_when(
      Compartment %in% c("Leaf", "Leaf endosphere") ~ "Leaf endosphere",
      Compartment %in% c("Root", "Root endosphere") ~ "Root endosphere",
      Compartment %in% c("Soil", "Rhizosphere")     ~ "Rhizosphere",
      TRUE                                          ~ as.character(Compartment))
  ) |>
  dplyr::filter(Compartment %in% c("Leaf endosphere", "Root endosphere", "Rhizosphere")) |>
  top_n_pathways(group_var = "Compartment", n = 15) |>
  dplyr::mutate(Compartment = factor(Compartment,
                                     levels = c("Leaf endosphere",
                                                "Root endosphere",
                                                "Rhizosphere")))

top_c <- top_compartment

# -------------------------------------------------------------------
# 7) Top 15 by drought (6 moisture classes)
# -------------------------------------------------------------------
# Create drought classes M1–M6 from watering amounts
pathway_drought <- pathway_cleaned |>
  dplyr::mutate(Drought_num = suppressWarnings(as.numeric(Drought)),
                Drought_class = dplyr::case_when(
                  Drought_num == 0                        ~ "M1: 0 ml",
                  Drought_num %in% c(300, 400)           ~ "M2: 300 & 400 ml",
                  Drought_num %in% c(600, 700)           ~ "M3: 600 & 700 ml",
                  Drought_num %in% c(900, 1000)          ~ "M4: 900 & 1000 ml",
                  Drought_num %in% c(1200, 1300)         ~ "M5: 1200 & 1300 ml",
                  Drought_num %in% c(1500, 1600)         ~ "M6: 1500 & 1600 ml",
                  TRUE                                   ~ NA_character_)) |>
  dplyr::filter(!is.na(Drought_class))

top_drought <- top_n_pathways(pathway_drought, group_var = "Drought_class", n = 15) |>
  dplyr::mutate(Drought_class = factor(Drought_class, levels = c("M1: 0 ml",
                                                                 "M2: 300 & 400 ml",
                                                                 "M3: 600 & 700 ml",
                                                                 "M4: 900 & 1000 ml",
                                                                 "M5: 1200 & 1300 ml",
                                                                 "M6: 1500 & 1600 ml")))

top_all <- top_drought |>
  dplyr::rename(Drought = Drought_class)

# -------------------------------------------------------------------
# 8) Ordering for level2 (KEGG level 2 category)
# -------------------------------------------------------------------
order_level2 <- c("Global and overview maps",
                  "Membrane transport",
                  "Signal transduction",
                  "Cellular community - prokaryotes",
                  "Nucleotide metabolism",
                  "Carbohydrate metabolism",
                  "Lipid metabolism")

top_all$level2 <- factor(top_all$level2, levels = order_level2)
top_c$level2   <- factor(top_c$level2,   levels = order_level2)
top_s$level2   <- factor(top_s$level2,   levels = order_level2)

# -------------------------------------------------------------------
# 9) Plotting
# -------------------------------------------------------------------
# Local theme (we keep it simple & consistent)
theme_pathway <- theme_classic() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 10, face = "bold"),
        axis.text.x  = element_text(size = 8,  face = "bold"),
        axis.text.y  = element_text(size = 10, face = "bold", color = "black"),
        axis.title   = element_text(size = 12, face = "bold"),
        strip.text   = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "#f0f4f8", color = NA))

# color palette for level2
col1 <- RColorBrewer::brewer.pal(8, "Set2")[1:7]

## 9.1 Drought (6 classes)
p_all <- ggplot(top_all, aes(x    = reorder(level1, Average_Abundance),
                             y    = Average_Abundance,
                             fill = level2)) +
  geom_col() +
  coord_flip() +
  theme_pathway +
  theme(strip.background = element_rect(fill = "#f3f9fe", color = NA)) +
  scale_fill_manual(values = col1) +
  facet_wrap(~ Drought, nrow = 2, ncol = 3) +
  labs(x    = "Level 3 KEGG pathway",
       y    = "Relative abundance (%)",
       fill = "Level 2 KEGG pathway")

## 9.2 Compartment
p_c <- ggplot(top_c, aes(x    = reorder(level1, Average_Abundance),
                         y    = Average_Abundance,
                         fill = level2)) +
  geom_col() +
  coord_flip() +
  theme_pathway +
  theme(strip.background = element_rect(fill = "#f9fef3", color = NA)) +
  scale_fill_manual(values = col1) +
  facet_wrap(~ Compartment, nrow = 1) +
  labs(x    = "Level 3 KEGG pathway",
       y    = "Relative abundance (%)",
       fill = "Level 2 KEGG pathway")

## 9.3 Plant species
p_s <- ggplot(top_s, aes(x    = reorder(level1, Average_Abundance),
                         y    = Average_Abundance,
                         fill = level2)) +
  geom_col() +
  coord_flip() +
  theme_pathway +
  theme(strip.text      = element_text(size = 10, face = "bold.italic"),
        legend.position = "right") +
  scale_fill_manual(values = col1) +
  facet_wrap(~ Species, nrow = 1) +
  labs(x    = "Level 3 KEGG pathway",
       y    = "Relative abundance (%)",
       fill = "Level 2 KEGG pathway")

# -------------------------------------------------------------------
# 10) Combined figure & export
# -------------------------------------------------------------------
f1 <- ggpubr::ggarrange(p_s, p_all, p_c, ncol = 1, heights = c(1, 2, 1),
                        labels = c("(a)", "(b)", "(c)"), font.label = list(size = 11))

ggsave(file.path(results_dir, "Top15_pathways_Bacteria.png"), ## Figure S12
       f1, width = 12, height = 12, units = "in", dpi = 600, bg = "white")

