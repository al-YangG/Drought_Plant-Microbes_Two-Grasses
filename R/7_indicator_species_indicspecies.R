##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/7_indicator_species_indicspecies.R
##
## Tasks:
##  1) Run indicspecies::multipatt for Bacteria & Fungi:
##       - Group_type = "Compartment", "Plant species", "Climatic conditions"
##  2) Extract significant indicators (p <= 0.001) and attach taxonomy
##  3) Build a combined indicator table for later use (DESeq2 overlap, etc.)
##  4) Summarise number of indicator ASVs per Group_type × Group × Domain
##  5) Build phylum-level indicator tables (Species / Compartment / Climate)
##  6) Plot:
##       (a) barplot of indicator counts per Group
##       (b–d) phylum-level bubble plots (Domain × phylum × Group)
##       (e) combined figure saved as results/Indicator_species.png
##
## Note:
##  - Uses phyloseq objects stored under Data/7_Indicator species
##  - Assumes taxonomy has phylum for both ps-Bac and ps-Fun_withGuilds
##  - Fungal phylum "Fungi_phy_Incertae_sedis" is recoded as "Unassigned"
## -------------------------------------------------------------------

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")
# (optional) if you later add helper functions:
# source("R/utils_indicator_taxa.R")

load_project_packages()

root      <- here::here()
data_dir  <- file.path(root, "Data", "7_Indicator species")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Load phyloseq objects & construct ASV + env tables
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(data_dir, "ps-Bac.rds"))
ps_Fun <- readRDS(file.path(data_dir, "ps-Fun_withGuilds.rds"))

# helper: convert phyloseq OTU table to matrix with samples in rows
otu_to_sample_matrix <- function(ps) {
  m <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

asv_B <- otu_to_sample_matrix(ps_Bac)  # Bacteria
asv_F <- otu_to_sample_matrix(ps_Fun)  # Fungi

env_B <- as.data.frame(phyloseq::sample_data(ps_Bac))
env_F <- as.data.frame(phyloseq::sample_data(ps_Fun))

# make sure key factors are correctly typed
for (df in list(env_B, env_F)) {
  df$Compartment <- factor(df$Compartment,
                           levels = c("Leaf endosphere",
                                      "Root endosphere",
                                      "Rhizosphere"))
  df$Species     <- factor(df$Species,
                           levels = c("Festuca", "Lolium"))
  df$Climate     <- factor(df$Climate,
                           levels = c("Current", "Future"))
}

env_B$Mesocosm <- factor(env_B$Mesocosm)
env_F$Mesocosm <- factor(env_F$Mesocosm)

# -------------------------------------------------------------------
# 2) Helper: run multipatt & extract significant indicators
# -------------------------------------------------------------------
# We treat three groupings:
#   - Group_type = "Compartment"         (env$Compartment)
#   - Group_type = "Plant species"       (env$Species)
#   - Group_type = "Climatic conditions" (env$Climate)
#
# For compartments we often exclude site-group combinations (duleg = TRUE).

run_indicators_one <- function(asv_mat,
                               groups,
                               group_type,
                               domain,
                               use_duleg = FALSE,
                               alpha_cut = 0.001,
                               nperm = 999) {
  stopifnot(nrow(asv_mat) == length(groups))
  
  # permutation design: use simple how() here
  ctrl <- permute::how(nperm = nperm)
  
  if (use_duleg) {
    fit <- indicspecies::multipatt(asv_mat, groups,
                                   duleg = TRUE, control = ctrl)
  } else {
    fit <- indicspecies::multipatt(asv_mat, groups,
                                   control = ctrl)
  }
  
  sign_tab <- as.data.frame(fit$sign)
  sign_tab$ASVid <- rownames(sign_tab)
  
  # keep only significant ASVs at alpha_cut
  sig <- subset(sign_tab, p.value <= alpha_cut)
  if (nrow(sig) == 0L) return(NULL)
  
  # For simple groups (not combinations) the columns named by group
  # contain 1/0 for membership. Identify the "group" for each ASV.
  group_cols <- colnames(sig)[colnames(sig) %in% levels(groups)]
  if (length(group_cols) == 0L) {
    stop("No group columns found in indicspecies output for ", group_type)
  }
  
  long <- sig |>
    tidyr::pivot_longer(cols = all_of(group_cols),
                        names_to = "Group", values_to = "Member") |>
    dplyr::filter(Member == 1) |>
    dplyr::select(ASVid, Group, stat = stat, p.value)
  
  if (!nrow(long)) return(NULL)
  
  long |>
    dplyr::mutate(Domain    = domain,
                  Group_type = group_type,
                  Group     = as.character(Group)) |>
    dplyr::relocate(Domain, Group_type, Group, ASVid)
}

# wrapper: run for all three groupings in one domain
run_indicators_domain <- function(asv_mat, env_df, domain_label) {
  dplyr::bind_rows(
    # Compartment (duleg = TRUE)
    run_indicators_one(asv_mat  = asv_mat,
                       groups   = env_df$Compartment,
                       group_type = "Compartment",
                       domain   = domain_label,
                       use_duleg = TRUE),
    # Plant species
    run_indicators_one(asv_mat  = asv_mat,
                       groups   = env_df$Species,
                       group_type = "Plant species",
                       domain   = domain_label,
                       use_duleg = FALSE),
    # Climatic conditions
    run_indicators_one(asv_mat  = asv_mat,
                       groups   = env_df$Climate,
                       group_type = "Climatic conditions",
                       domain   = domain_label,
                       use_duleg = FALSE))
}

# -------------------------------------------------------------------
# 3) Run indicspecies for Bacteria and Fungi, combine results
# -------------------------------------------------------------------
ind_B <- run_indicators_domain(asv_B, env_B, domain_label = "Bacteria")
ind_F <- run_indicators_domain(asv_F, env_F, domain_label = "Fungi")

ind_all_raw <- dplyr::bind_rows(ind_B, ind_F)

# if nothing is significant, bail early
if (is.null(ind_all_raw) || !nrow(ind_all_raw)) {
  stop("No significant indicators found at p <= 0.001 – check inputs.")
}

# -------------------------------------------------------------------
# 4) Attach taxonomy (to phylum level) and clean fungal phylum names
# -------------------------------------------------------------------
tax_B <- as.data.frame(phyloseq::tax_table(ps_Bac)) |>
  tibble::rownames_to_column("ASVid")

tax_F <- as.data.frame(phyloseq::tax_table(ps_Fun)) |>
  tibble::rownames_to_column("ASVid") |>
  dplyr::mutate(phylum = dplyr::if_else(phylum == "Fungi_phy_Incertae_sedis",
                                        "Unassigned",
                                        as.character(phylum)))

ind_with_tax <- ind_all_raw |>
  dplyr::left_join(tax_B, by = "ASVid") |>
  dplyr::mutate(dplyr::across(
    dplyr::all_of(colnames(tax_F)[colnames(tax_F) != "ASVid"]),
    ~ dplyr::coalesce(.x, tax_F[[cur_column()]][match(ASVid, tax_F$ASVid)])))

# ensure phylum column exists and is character
ind_with_tax$phylum <- as.character(ind_with_tax$phylum)

# tidy names / columns
indicspecies_tbl <- ind_with_tax |>
  dplyr::rename(Stat = stat, P_value = p.value) |>
  dplyr::arrange(Domain, Group_type, Group, P_value)

# save combined table for later (DESeq2 overlap etc.)
write.csv(indicspecies_tbl,
          file.path(results_dir, "All_indicator_species_indicspecies.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 5) Summary: number of indicator ASVs per Group_type × Group × Domain
# -------------------------------------------------------------------
indicator_number_tbl <- indicspecies_tbl |>
  dplyr::count(Group_type, Group, Domain, name = "Number") |>
  dplyr::arrange(Group_type, Group, Domain)

write.csv(indicator_number_tbl,
          file.path(results_dir, "Indicator_number.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 6) Build phylum-level indicator tables (for bubble plots)
# -------------------------------------------------------------------
# total ASVs per Domain × phylum (for denominator)
tot_bac <- tax_B |>
  dplyr::count(phylum, name = "Total_ASVs") |>
  dplyr::mutate(Domain = "Bacteria")

tot_fun <- tax_F |>
  dplyr::count(phylum, name = "Total_ASVs") |>
  dplyr::mutate(Domain = "Fungi")

tot_all <- dplyr::bind_rows(tot_bac, tot_fun)

# helper to summarise indicators per Group_type × Group × phylum
summarise_phylum_indicators <- function(df, group_type_name) {
  df |>
    dplyr::filter(Group_type == group_type_name) |>
    dplyr::group_by(Domain, {{ group_type_name := Group }}, phylum) |>
    dplyr::summarise(n_ASVs   = dplyr::n(),
                     meanStat = mean(Stat, na.rm = TRUE),
                     .groups  = "drop")
}

# For clarity we do it explicitly:
ind_comp <- indicspecies_tbl |>
  dplyr::filter(Group_type == "Compartment") |>
  dplyr::group_by(Domain, Compartment = Group, phylum) |>
  dplyr::summarise(n_ASVs   = dplyr::n(),
                   meanStat = mean(Stat, na.rm = TRUE),
                   .groups  = "drop")

ind_sp <- indicspecies_tbl |>
  dplyr::filter(Group_type == "Plant species") |>
  dplyr::group_by(Domain, Species = Group, phylum) |>
  dplyr::summarise(n_ASVs   = dplyr::n(),
                   meanStat = mean(Stat, na.rm = TRUE),
                   .groups  = "drop")

ind_clim <- indicspecies_tbl |>
  dplyr::filter(Group_type == "Climatic conditions") |>
  dplyr::group_by(Domain, Climate = Group, phylum) |>
  dplyr::summarise(n_ASVs   = dplyr::n(),
                   meanStat = mean(Stat, na.rm = TRUE),
                   .groups  = "drop")

# add total indicator counts and total ASVs for labels

## species
phylum_tot_sp <- ind_sp |>
  dplyr::group_by(Domain, phylum) |>
  dplyr::summarise(n_ASVs_total = sum(n_ASVs),
                   .groups = "drop")

ind_sp_lab <- ind_sp |>
  dplyr::left_join(phylum_tot_sp, by = c("Domain", "phylum")) |>
  dplyr::left_join(tot_all,       by = c("Domain", "phylum")) |>
  dplyr::mutate(phylum_label = paste0(phylum, "\n(", n_ASVs_total, " / ", Total_ASVs, ")"))

## compartments
phylum_tot_comp <- ind_comp |>
  dplyr::group_by(Domain, phylum) |>
  dplyr::summarise(total_ind = sum(n_ASVs), .groups   = "drop")

ind_comp_lab <- ind_comp |>
  dplyr::left_join(phylum_tot_comp, by = c("Domain", "phylum")) |>
  dplyr::left_join(tot_all,         by = c("Domain", "phylum")) |>
  dplyr::mutate(phylum_label = paste0(phylum, "\n(", total_ind, " / ", Total_ASVs, ")"))

## climates
phylum_tot_clim <- ind_clim |>
  dplyr::group_by(Domain, phylum) |>
  dplyr::summarise(total_ind = sum(n_ASVs), .groups   = "drop")

ind_clim_lab <- ind_clim |>
  dplyr::left_join(phylum_tot_clim, by = c("Domain", "phylum")) |>
  dplyr::left_join(tot_all,         by = c("Domain", "phylum")) |>
  dplyr::mutate(phylum_label = paste0(phylum, "\n(", total_ind, " / ", Total_ASVs, ")"))

# order factors
ind_comp_lab$Compartment <- factor(
  ind_comp_lab$Compartment,
  levels = c("Leaf endosphere", "Root endosphere", "Rhizosphere")
)

ind_clim_lab$Climate <- factor(
  ind_clim_lab$Climate,
  levels = c("Current", "Future")
)

ind_sp_lab$Species <- factor(
  ind_sp_lab$Species,
  levels = c("Festuca", "Lolium")
)

# reorder phylum labels by total indicator counts within domain
ind_sp_lab <- ind_sp_lab |>
  dplyr::group_by(phylum_label) |>
  dplyr::mutate(total_phylum = dplyr::first(n_ASVs_total)) |>
  dplyr::ungroup() |>
  dplyr::mutate(phylum_label = reorder(phylum_label, total_phylum))

ind_comp_lab <- ind_comp_lab |>
  dplyr::mutate(total_phylum = total_ind,
                phylum_label = reorder(phylum_label, total_phylum))

ind_clim_lab <- ind_clim_lab |>
  dplyr::mutate(total_phylum = total_ind,
                phylum_label = reorder(phylum_label, total_phylum))

# optionally save the phylum-level indicator tables
write.csv(ind_sp_lab,
          file.path(results_dir, "Indicator_taxa_bySpecies.csv"),
          row.names = FALSE)

write.csv(ind_comp_lab,
          file.path(results_dir, "Indicator_taxa_byCompartment.csv"),
          row.names = FALSE)

write.csv(ind_clim_lab,
          file.path(results_dir, "Indicator_taxa_byClimate.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 7) Plotting – overview numbers & phylum-level bubble plots
# -------------------------------------------------------------------

# --- overview bar plot: indicator counts per Group_type × Group × Domain
indicator_number_plot <- {
  ind_clean <- indicator_number_tbl |>
    dplyr::mutate(Group_type = factor(
      Group_type, levels = c("Compartment", "Plant species", "Climatic conditions")),
      Domain = factor(Domain, levels = c("Bacteria", "Fungi")),
      Group  = factor(Group, levels = c("Leaf endosphere",
                                        "Root endosphere",
                                        "Rhizosphere",
                                        "Festuca",
                                        "Lolium",
                                        "Current",
                                        "Future")))
  
  dodge <- position_dodge(width = 0.7)
  
  ggplot(ind_clean, aes(x = Group, y = Number, fill = Domain)) +
    geom_col(position = dodge, width = 0.6) +
    geom_text(aes(label = Number), position = dodge, vjust = -0.5, size = 3.5, fontface = "bold") +
    facet_wrap(~ Group_type, scales = "free_x") +
    scale_fill_manual(values = c("Bacteria" = "#1B9E77", "Fungi" = "#A6761D")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of indicator ASVs", fill = "Domain") +
    theme_classic() +
    theme(strip.background = element_rect(fill = "grey95", color = "white"),
          strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
          axis.title.y = element_text(face = "bold", size = 11),
          panel.grid.major.x = element_blank(),
          panel.spacing = grid::unit(1.2, "lines"),
          legend.justification = c(0, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.text = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 11))
}

# theme for bubble plots
indicator_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(face = "bold", size = 11, angle = 30, hjust = 1, color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10),
        panel.grid.major.x = element_blank())

# --- bubble plot: Species × phylum
p_species <- ggplot(ind_sp_lab, aes(x = Species, y = phylum_label)) +
  geom_point(aes(size = n_ASVs, color = meanStat)) +
  facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Number of\nindicator ASVs") +
  scale_color_viridis_c(option = "plasma", direction = -1,name = "Mean IndVal") +
  labs(x = NULL, y = "Phylum (indicator / total ASVs)") +
  indicator_theme +
  theme(axis.text.x = element_text(face = "bold.italic", size = 11, color = "black"),
        legend.position = "right") +
  guides(size = guide_legend(order = 1))

# --- bubble plot: Compartment × phylum
p_comp <- ggplot(ind_comp_lab, aes(x = Compartment, y = phylum_label)) +
  geom_point(aes(size = n_ASVs, color = meanStat)) +
  facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Number of\nindicator ASVs") +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Mean IndVal") +
  labs(x = NULL, y = "Phylum (indicator / total ASVs)") +
  indicator_theme +
  theme(axis.text.x = element_text(face = "bold", size = 11, angle = 30, hjust = 1, color = "black"),
        legend.position = "right") +
  guides(size = guide_legend(order = 1))

# --- bubble plot: Climate × phylum
p_clim <- ggplot(ind_clim_lab, aes(x = Climate, y = phylum_label)) +
  geom_point(aes(size = n_ASVs, color = meanStat)) +
  facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Number of\nindicator ASVs") +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Mean IndVal") +
  labs(x = NULL, y = "Phylum (indicator / total ASVs)") +
  indicator_theme +
  theme(axis.text.x = element_text(face = "bold", size = 11, hjust = 0.5, color = "black"),
        legend.position = "right") +
  guides(size = guide_legend(order = 1))

# -------------------------------------------------------------------
# 8) Assemble composite figure and save
# -------------------------------------------------------------------
p_ind   <- indicator_number_plot
p_three <- ggpubr::ggarrange(p_comp, p_species, p_clim, nrow = 1,
                             labels = c("(b)", "(c)", "(d)"),
                             widths = c(1.15, 1, 1))

final_fig <- ggpubr::ggarrange(p_ind, p_three, ncol = 1,
                               labels = c("(a)", " "), heights = c(1, 2))

ggsave(filename = file.path(results_dir, "Indicator_species.png"), ## Figure S16
       plot = final_fig, width = 14, height = 10, units = "in", dpi = 600, bg = "white")

