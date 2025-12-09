##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/7c_indicator_DESeq2_summaries_vs_indicspecies.R
##
## Tasks:
##  1) Run DESeq2 moisture-gradient models for:
##       - Plant species (Festuca, Lolium)
##       - Compartment
##       - Climatic conditions
##     separately for Bacteria (ps-Bac) and Fungi (ps-Fun_withGuilds).
##  2) Summarise counts of moisture-responsive ASVs (Positive/Negative)
##     per Domain × Group_type × Group.
##  3) Visualise DESeq2 counts as stacked bars.
##  4) Compare indicspecies vs DESeq2:
##       - Overlap counts (Indicspecies only / Both / DESeq2 only)
##       - Summary plot across Compartment, Plant species, Climate.
##  5) Export overlap table and common ASVs with basic taxonomy.
##
## Note:
##  - This script *re-runs* DESeq2 grouped analyses. If you prefer to reuse
##    results from 7b, you can comment out the DESeq section and just
##    read pre-computed CSVs.
## -------------------------------------------------------------------

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")

load_project_packages()

root        <- here::here()
data_dir    <- file.path(root, "Data", "7_Indicator species")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Load phyloseq objects & taxonomy
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(data_dir, "ps-Bac.rds"))
ps_Fun <- readRDS(file.path(data_dir, "ps-Fun_withGuilds.rds"))

# taxonomy as data.frame with ASV
tax_Bac <- phyloseq::tax_table(ps_Bac) |>
  as.data.frame() |>
  tibble::rownames_to_column("ASV")

tax_Fun <- phyloseq::tax_table(ps_Fun) |>
  as.data.frame() |>
  tibble::rownames_to_column("ASV") |>
  dplyr::mutate(phylum = dplyr::if_else(phylum == "Fungi_phy_Incertae_sedis",
                                        "Unassigned", as.character(phylum)))

# ensure scaled Moisture column exists in both ps objects
for (ps in list(ps_Bac, ps_Fun)) {
  sd <- phyloseq::sample_data(ps)
  if (!"Moisture_sc" %in% colnames(sd)) {
    sd$Moisture_sc <- as.numeric(scale(sd$Moisture))
  } else {
    sd$Moisture_sc <- as.numeric(sd$Moisture_sc)
  }
  phyloseq::sample_data(ps) <- sd
}
ps_Bac <- ps_Bac
ps_Fun <- ps_Fun

# -------------------------------------------------------------------
# 2) Generic DESeq2 helpers
# -------------------------------------------------------------------
# classify gradient sign & significance
classify_trend <- function(res_df, p_cut = 0.05) {
  res_df$Diff.abund <- "Not sig"
  res_df$Diff.abund[res_df$pvalue < p_cut & res_df$log2FoldChange > 0] <- "Positive"
  res_df$Diff.abund[res_df$pvalue < p_cut & res_df$log2FoldChange < 0] <- "Negative"
  res_df
}

# Generic helper:
#  run DESeq2 gradient within levels of a grouping variable
run_deseq_gradient_group <- function(ps, tax_df, domain_label,
                                     group_var,      # e.g. "Species", "Compartment", "Climate"
                                     group_levels,   # character vector of levels to run
                                     group_type_lab  # e.g. "Plant species", "Compartment"
) {
  purrr::map_dfr(group_levels, function(g_lvl) {
    message("Running DESeq2 for ", domain_label,
            " | ", group_type_lab, " = ", g_lvl)
    
    sd <- phyloseq::sample_data(ps)
    if (!group_var %in% colnames(sd)) {
      warning("Grouping variable ", group_var, " not found in sample_data.")
      return(NULL)
    }
    
    keep <- sd[[group_var]] == g_lvl
    keep[is.na(keep)] <- FALSE
    
    ps_sub <- phyloseq::prune_samples(keep, ps)
    if (phyloseq::nsamples(ps_sub) == 0L) {
      warning("No samples for ", group_var, " = ", g_lvl,
              " in ", domain_label)
      return(NULL)
    }
    
    # ensure Moisture_sc
    sd_sub <- phyloseq::sample_data(ps_sub)
    if (!"Moisture_sc" %in% colnames(sd_sub)) {
      sd_sub$Moisture_sc <- as.numeric(scale(sd_sub$Moisture))
      phyloseq::sample_data(ps_sub) <- sd_sub
    } else {
      sd_sub$Moisture_sc <- as.numeric(sd_sub$Moisture_sc)
      phyloseq::sample_data(ps_sub) <- sd_sub
    }
    
    # build DESeq2 object
    dds <- phyloseq::phyloseq_to_deseq2(ps_sub, ~ Moisture_sc)
    
    # drop zero libraries
    lib_sizes <- colSums(DESeq2::counts(dds))
    keep_lib  <- lib_sizes > 0
    dds       <- dds[, keep_lib]
    
    # normalise
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric")
    
    # extract moisture gradient effect
    res <- as.data.frame(DESeq2::results(dds, cooksCutoff = FALSE))
    res <- classify_trend(res)
    res$ASV <- rownames(res)
    
    res_tax <- dplyr::left_join(res, tax_df, by = "ASV")
    
    res_tax <- res_tax |>
      dplyr::mutate(Domain = domain_label, Group_type = group_type_lab, Group = g_lvl)
    
    res_tax
  })
}

# -------------------------------------------------------------------
# 3) DESeq2: species, compartment, climate (Bacteria + Fungi)
# -------------------------------------------------------------------

## 3.1 plant species
species_levels <- sort(unique(as.character(phyloseq::sample_data(ps_Bac)$Species)))
species_levels <- species_levels[!is.na(species_levels)]

deseq_sp_Bac <- run_deseq_gradient_group(ps = ps_Bac, tax_df = tax_Bac,
                                         domain_label = "Bacteria", group_var = "Species",
                                         group_levels = species_levels,
                                         group_type_lab = "Plant species")

deseq_sp_Fun <- run_deseq_gradient_group(ps = ps_Fun, tax_df = tax_Fun,
                                         domain_label = "Fungi", group_var = "Species",
                                         group_levels = species_levels, 
                                         group_type_lab = "Plant species")

deseq_sp_all <- dplyr::bind_rows(deseq_sp_Bac, deseq_sp_Fun)
deseq_sp_sig <- dplyr::filter(deseq_sp_all, Diff.abund != "Not sig")

write.csv(deseq_sp_sig, row.names = FALSE,
          file.path(results_dir, "DESeq_Bac_Fun_Moisture_bySpecies.csv"))

## 3.2 compartment (use whatever levels exist in ps_Bac)
comp_levels <- sort(unique(as.character(phyloseq::sample_data(ps_Bac)$Compartment)))
comp_levels <- comp_levels[!is.na(comp_levels)]

deseq_comp_Bac <- run_deseq_gradient_group(ps = ps_Bac, tax_df = tax_Bac,
                                           domain_label = "Bacteria",
                                           group_var = "Compartment",
                                           group_levels = comp_levels,
                                           group_type_lab = "Compartment")

deseq_comp_Fun <- run_deseq_gradient_group(ps = ps_Fun, tax_df = tax_Fun,
                                           domain_label = "Fungi",
                                           group_var = "Compartment",
                                           group_levels = comp_levels,
                                           group_type_lab = "Compartment")

deseq_comp_all <- dplyr::bind_rows(deseq_comp_Bac, deseq_comp_Fun)
deseq_comp_sig <- dplyr::filter(deseq_comp_all, Diff.abund != "Not sig")

write.csv(deseq_comp_sig, row.names = FALSE,
          file.path(results_dir, "DESeq_Bac_Fun_Moisture_byCompartment.csv"))

## 3.3 climate
clim_levels <- sort(unique(as.character(phyloseq::sample_data(ps_Bac)$Climate)))
clim_levels <- clim_levels[!is.na(clim_levels)]

deseq_clim_Bac <- run_deseq_gradient_group(ps = ps_Bac, tax_df = tax_Bac,
                                           domain_label = "Bacteria",
                                           group_var = "Climate",
                                           group_levels = clim_levels,
                                           group_type_lab = "Climatic conditions")

deseq_clim_Fun <- run_deseq_gradient_group(ps = ps_Fun, tax_df = tax_Fun,
                                           domain_label = "Fungi",
                                           group_var = "Climate",
                                           group_levels = clim_levels,
                                           group_type_lab = "Climatic conditions")

deseq_clim_all <- dplyr::bind_rows(deseq_clim_Bac, deseq_clim_Fun)
deseq_clim_sig <- dplyr::filter(deseq_clim_all, Diff.abund != "Not sig")

write.csv(deseq_clim_sig, row.names = FALSE,
          file.path(results_dir, "DESeq_Bac_Fun_Moisture_byClimate.csv"))

# -------------------------------------------------------------------
# 4) Summary: number of DESeq2-responsive ASVs per group
# -------------------------------------------------------------------
# species-level
summ_species <- deseq_sp_sig |>
  dplyr::count(Domain, Group_type, Group, Diff.abund, name = "n")

# compartment-level
summ_compartment <- deseq_comp_sig |>
  dplyr::count(Domain, Group_type, Group, Diff.abund, name = "n")

# climate-level
summ_climate <- deseq_clim_sig |>
  dplyr::count(Domain, Group_type, Group, Diff.abund, name = "n")

summ_all <- dplyr::bind_rows(summ_species, summ_compartment, summ_climate)

write.csv(summ_all, row.names = FALSE,
          file.path(results_dir, "DESeq_number_per_group.csv"))

# nice factor ordering for plotting
summ_all <- summ_all |>
  dplyr::mutate(Group_type = factor(Group_type,
                                    levels = c("Compartment",
                                               "Plant species",
                                               "Climatic conditions")),
    Domain = factor(Domain, levels = c("Bacteria", "Fungi")),
    Diff.abund = factor(Diff.abund, levels = c("Positive", "Negative")),
    Group = factor(Group, levels = c("Leaf endosphere",
                                     "Root endosphere",
                                     "Rhizosphere",
                                     "Festuca",
                                     "Lolium",
                                     "Current",
                                     "Future")))

p_summary <- ggplot(summ_all, aes(x = Group, y = n, fill = Diff.abund)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  facet_grid(Domain ~ Group_type, scales = "free_x", space  = "free_x") +
  labs(x = NULL, y = "Number of moisture-responsive ASVs (DESeq2)", fill = "Direction") +
  scale_fill_manual(values = c("Positive" = "#fb9a99","Negative" = "#a6cee3")) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = -30, hjust = 0, face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        legend.position = c(0.80, 0.88),
        legend.justification = c(0, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11))

# -------------------------------------------------------------------
# 5) indicspecies vs DESeq2 overlap
# -------------------------------------------------------------------
df_ind <- read.csv(file.path(data_dir, "All_Indicator species_indicspecies.csv"))
df_DES <- read.csv(file.path(data_dir, "All_Indicator species_DESeq2.csv"))

df_ind <- df_ind |>
  dplyr::mutate(ASVid = as.character(ASVid))
df_DES <- df_DES |>
  dplyr::mutate(ASVid = as.character(ASVid))

# helper: ASV sets per combination
get_sets <- function(gt, g, d, df_ind, df_DES) {
  ind_ids <- df_ind |>
    dplyr::filter(Group_type == gt, Group == g, Domain == d) |>
    dplyr::pull(ASVid) |>
    unique()
  
  des_ids <- df_DES |>
    dplyr::filter(Group_type == gt, Group == g, Domain == d) |>
    dplyr::pull(ASVid) |>
    unique()
  
  list(ind = ind_ids, des = des_ids)
}

# helper: overlap counts
sum_overlap <- function(sets) {
  ind <- sets$ind
  des <- sets$des
  
  tibble::tibble(ind_only = sum(!ind %in% des),
                 des_only = sum(!des %in% ind),
                 both     = sum(ind %in% des))
}

# all unique combinations
combs <- df_ind |>
  dplyr::select(Group_type, Group, Domain) |>
  dplyr::distinct()

overlap_tbl <- combs |>
  dplyr::mutate(sets   = purrr::pmap(list(Group_type, Group, Domain),
                                     ~ get_sets(..1, ..2, ..3, df_ind, df_DES)),
                counts = purrr::map(sets, sum_overlap)) |>
  dplyr::select(-sets) |>
  tidyr::unnest(counts) |>
  tidyr::pivot_longer(cols      = c(ind_only, des_only, both),
                      names_to  = "Category",
                      values_to = "Count") |>
  dplyr::mutate(Category = factor(Category,
                                  levels = c("ind_only", "both", "des_only"),
                                  labels = c("Indicspecies only", "Both", "DESeq2 only")))

# separate tables by Group_type to control ordering if needed
comp_order <- c("Leaf endosphere", "Root endosphere", "Rhizosphere")

dat_clim <- overlap_tbl |>
  dplyr::filter(Group_type == "Climatic conditions")

dat_comp <- overlap_tbl |>
  dplyr::filter(Group_type == "Compartment") |>
  dplyr::mutate(Group = factor(Group, levels = comp_order))

dat_sp <- overlap_tbl |>
  dplyr::filter(Group_type == "Plant species")

plot_overlap_gt <- function(data) {
  data$Group <- factor(data$Group)
  pos <- position_dodge(width = 0.6)
  
  ggplot(data, aes(x = Group, y = Count, fill = Category)) +
    geom_point(size = 8, shape = 21, color = "black", stroke = 0.4, position = pos) +
    geom_text(aes(label = Count), color = "white", fontface = "bold", size = 3,
              position = pos, show.legend = FALSE) +
    facet_grid(. ~ Domain, scales = "free_y", space  = "free_y") +
    scale_y_continuous(trans  = "sqrt", limits = c(0, NA),
                       breaks = c(0, 10, 50, 100, 250, 500, 1000),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("Indicspecies only" = "#3288BD",
                                 "Both"              = "#E41A1C",
                                 "DESeq2 only"       = "#5E4FA2")) +
    labs(x = NULL, y = "Number of ASVs", fill = "Category") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, color = "black", face = "bold"),
          axis.text.y = element_text(color = "black", face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold", size = 11))
}

p_clim <- plot_overlap_gt(dat_clim)
p_comp <- plot_overlap_gt(dat_comp)
p_sp   <- plot_overlap_gt(dat_sp)

p_num <- ggpubr::ggarrange(p_comp, p_sp, p_clim, ncol = 1, heights = c(1.3, 1, 1),
                           labels = c("(b)", "(c)", "(d)"), common.legend = TRUE, legend = "bottom")

p1 <- ggpubr::ggarrange(p_summary, p_num, nrow = 1, labels = c("(a)", ""))

ggsave(file.path(results_dir, "Indicator_species_indicspecies_vs_DESeq2.png"), ## Figure S14
       p1, width = 14, height = 8, units = "in", dpi = 600, bg = "white")

# -------------------------------------------------------------------
# 6) Detailed common ASVs (indicspecies ∩ DESeq2)
# -------------------------------------------------------------------
get_common_details <- function(gt, g, d, df_ind, df_DES) {
  ind_sub <- df_ind |>
    dplyr::filter(Group_type == gt, Group == g, Domain == d)
  
  des_sub <- df_DES |>
    dplyr::filter(Group_type == gt, Group == g, Domain == d)
  
  common_ids <- intersect(unique(ind_sub$ASVid), unique(des_sub$ASVid))
  
  if (length(common_ids) == 0L) {
    return(tibble::tibble())
  }
  
  ind_common <- ind_sub |>
    dplyr::filter(ASVid %in% common_ids)
  
  des_common <- des_sub |>
    dplyr::filter(ASVid %in% common_ids)
  
  dplyr::full_join(ind_common, des_common,
                   by = c("ASVid", "Group_type", "Group", "Domain"),
                   suffix = c("_indicspecies", "_DESeq2")) |>
    dplyr::select(-Group_type, -Group, -Domain)
}

common_ASV_details <- combs |>
  dplyr::mutate(common = purrr::pmap(list(Group_type, Group, Domain),
                                     ~ get_common_details(..1, ..2, ..3, df_ind, df_DES))) |>
  tidyr::unnest(common)

common_ASV_details_clean <- common_ASV_details |>
  dplyr::mutate(phylum  = dplyr::coalesce(phylum_indicspecies,  phylum_DESeq2),
                class   = dplyr::coalesce(class_indicspecies,   class_DESeq2),
                order   = dplyr::coalesce(order_indicspecies,   order_DESeq2),
                family  = dplyr::coalesce(family_indicspecies,  family_DESeq2),
                genus   = dplyr::coalesce(genus_indicspecies,   genus_DESeq2),
                species = dplyr::coalesce(species_indicspecies, species_DESeq2),
                TrophicMode = dplyr::coalesce(TrophicMode_indicspecies, TrophicMode_DESeq2)) |>
  dplyr::select(Group_type, Group, Domain, ASVid, Stat, Direction, baseMean, log2FoldChange,
                lfcSE, stat, pvalue, phylum, class, order, family, genus, species, TrophicMode)

write.csv(common_ASV_details_clean, row.names = FALSE,
          file.path(results_dir, "Indicator_common_ASVs_indicspecies_DESeq2.csv")) ## Table S13

