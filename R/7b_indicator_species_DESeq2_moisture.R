##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/7b_indicator_species_DESeq2_moisture.R
##
## Tasks:
##  1) Run DESeq2 moisture-gradient models for:
##       - All Bacteria (ps-Bac)
##       - All Fungi    (ps-Fun_withGuilds)
##  2) Run moisture-gradient DESeq2 separately for:
##       - Festuca vs. Lolium, for both domains
##  3) Classify ASVs as Positive / Negative / Not sig along the moisture gradient
##  4) Save:
##       - Domain-level gradient results (Bacteria / Fungi)
##       - Species-specific moisture-responsive ASVs across domains
##  5) Visualise:
##       (a) Positive vs Negative counts per Species × Domain
##       (b) Species overlap (pseudo-venn) and shared ASV responses
##       (c) Lollipop plots by phylum for Bacteria and Fungi
##       (d) Combined figures saved to results/
## -------------------------------------------------------------------

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
source("R/utils_packages.R")
source("R/utils_theme.R")
# (optional) indicator helpers if you later add them:
# source("R/utils_indicator_taxa.R")

load_project_packages()

root        <- here::here()
data_dir    <- file.path(root, "Data", "7_Indicator species")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Load phyloseq objects & ensure scaled Moisture
# -------------------------------------------------------------------
ps_Bac <- readRDS(file.path(data_dir, "ps-Bac.rds"))
ps_Fun <- readRDS(file.path(data_dir, "ps-Fun_withGuilds.rds"))

# ensure Moisture_sc exists and is numeric
for (ps in list(ps_Bac, ps_Fun)) {
  sd <- phyloseq::sample_data(ps)
  if (!"Moisture_sc" %in% colnames(sd)) {
    sd$Moisture_sc <- as.numeric(scale(sd$Moisture))
  } else {
    sd$Moisture_sc <- as.numeric(sd$Moisture_sc)
  }
  phyloseq::sample_data(ps) <- sd
}
# reassign (list loop modifies copies)
ps_Bac <- ps_Bac
ps_Fun <- ps_Fun

# taxonomy data frames with ASV IDs
tax_Bac <- phyloseq::tax_table(ps_Bac) |>
  as.data.frame() |>
  tibble::rownames_to_column("ASV")

tax_Fun <- phyloseq::tax_table(ps_Fun) |>
  as.data.frame() |>
  tibble::rownames_to_column("ASV") |>
  dplyr::mutate(
    phylum = dplyr::if_else(
      phylum == "Fungi_phy_Incertae_sedis",
      "Unassigned",
      as.character(phylum)
    )
  )

# -------------------------------------------------------------------
# 2) Helper functions
# -------------------------------------------------------------------
# classify gradient sign & significance
classify_trend <- function(res_df, p_cut = 0.05) {
  res_df$Diff.abund <- "Not sig"
  res_df$Diff.abund[res_df$pvalue < p_cut & res_df$log2FoldChange > 0] <- "Positive"
  res_df$Diff.abund[res_df$pvalue < p_cut & res_df$log2FoldChange < 0] <- "Negative"
  res_df
}

# run DESeq2 moisture-gradient for full domain (all samples)
run_deseq_gradient_domain <- function(ps, tax_df, domain_label) {
  message("Running domain-level moisture gradient for ", domain_label)
  
  dds <- phyloseq::phyloseq_to_deseq2(ps, ~ Moisture_sc)
  
  # drop samples with zero library size
  lib_sizes <- colSums(DESeq2::counts(dds))
  keep      <- lib_sizes > 0
  dds       <- dds[, keep]
  
  # library-size normalisation
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  
  # DESeq
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric")
  
  # extract results
  res <- as.data.frame(DESeq2::results(dds, cooksCutoff = FALSE))
  res <- classify_trend(res)
  res$ASV <- rownames(res)
  
  res_tax <- dplyr::left_join(res, tax_df, by = "ASV")
  res_tax$Domain <- domain_label
  
  res_tax
}

# run DESeq2 moisture-gradient for one plant species
run_deseq_gradient_species <- function(ps, species_name, tax_df, domain_label) {
  message("Running species-level moisture gradient for ", domain_label,
          ", Species = ", species_name)
  
  sd <- phyloseq::sample_data(ps)
  keep_species <- sd$Species == species_name
  keep_species[is.na(keep_species)] <- FALSE
  
  ps_sub <- phyloseq::prune_samples(keep_species, ps)
  if (phyloseq::nsamples(ps_sub) == 0L) {
    warning("No samples for Species = ", species_name, " in ", domain_label)
    return(list(all = NULL, sig = NULL))
  }
  
  sd_sub <- phyloseq::sample_data(ps_sub)
  if (!"Moisture_sc" %in% colnames(sd_sub)) {
    sd_sub$Moisture_sc <- as.numeric(scale(sd_sub$Moisture))
    phyloseq::sample_data(ps_sub) <- sd_sub
  } else {
    sd_sub$Moisture_sc <- as.numeric(sd_sub$Moisture_sc)
    phyloseq::sample_data(ps_sub) <- sd_sub
  }
  
  dds <- phyloseq::phyloseq_to_deseq2(ps_sub, ~ Moisture_sc)
  
  lib_sizes <- colSums(DESeq2::counts(dds))
  keep      <- lib_sizes > 0
  dds       <- dds[, keep]
  
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric")
  
  res <- as.data.frame(DESeq2::results(dds, cooksCutoff = FALSE))
  res <- classify_trend(res)
  res$ASV <- rownames(res)
  
  res_tax <- dplyr::left_join(res, tax_df, by = "ASV")
  res_tax$Species_group <- species_name
  res_tax$Domain        <- domain_label
  
  sig <- res_tax[res_tax$Diff.abund != "Not sig", , drop = FALSE]
  sig <- sig[order(sig$pvalue), , drop = FALSE]
  
  message("  -> ", nrow(sig), " moisture-responsive ASVs")
  list(all = res_tax, sig = sig)
}

# -------------------------------------------------------------------
# 3) Domain-level DESeq (all samples) – Bacteria & Fungi
# -------------------------------------------------------------------
deseq_Bac_all <- run_deseq_gradient_domain(ps_Bac, tax_Bac, domain_label = "Bacteria")
deseq_Fun_all <- run_deseq_gradient_domain(ps_Fun, tax_Fun, domain_label = "Fungi")

# keep only significant moisture-responsive ASVs if desired
bac_sig_all_domain <- subset(deseq_Bac_all, Diff.abund != "Not sig")
fun_sig_all_domain <- subset(deseq_Fun_all, Diff.abund != "Not sig")

write.csv(bac_sig_all_domain, row.names = FALSE,
          file.path(results_dir, "DESeq_Bac_Moisture.csv"))

write.csv(fun_sig_all_domain, row.names = FALSE,
          file.path(results_dir, "DESeq_Fun_Moisture.csv"))

# -------------------------------------------------------------------
# 4) Species-level DESeq (Festuca & Lolium) – Bacteria & Fungi
# -------------------------------------------------------------------
# Bacteria
bac_Festuca <- run_deseq_gradient_species(ps_Bac, "Festuca", tax_Bac, "Bacteria")
bac_Lolium  <- run_deseq_gradient_species(ps_Bac, "Lolium",  tax_Bac, "Bacteria")

# Fungi
fun_Festuca <- run_deseq_gradient_species(ps_Fun, "Festuca", tax_Fun, "Fungi")
fun_Lolium  <- run_deseq_gradient_species(ps_Fun, "Lolium",  tax_Fun, "Fungi")

# bind significant ASVs from all species × domains
bac_sig_species <- dplyr::bind_rows(bac_Festuca$sig, bac_Lolium$sig)

fun_sig_species <- dplyr::bind_rows(fun_Festuca$sig, fun_Lolium$sig)

all_sig <- dplyr::bind_rows(bac_sig_species, fun_sig_species)

# simple sanity table
with(all_sig, table(Domain, Species_group, Diff.abund))

# save combined species-level DESeq table
write.csv(all_sig, row.names = FALSE,
          file.path(results_dir, "DESeq_Bac_Fun_Moisture_bySpecies.csv")) ## Table S12

# -------------------------------------------------------------------
# 5) Visualisation: direction barplot & species overlap (pseudo-venn)
# -------------------------------------------------------------------

# ---- 5.1 Positive vs Negative per Species × Domain ----
counts_dir <- all_sig |>
  dplyr::count(Domain, Species_group, Diff.abund)

theme_de_bar <- theme_classic() +
  theme(legend.position = c(0.95, 0.85),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 11, color = "black", face = "bold.italic"),
        axis.title = element_text(size = 11, face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(color = NA, fill = "gray90"))

p_dir_bar <- ggplot(counts_dir, aes(x = Species_group, y = n, fill = Diff.abund)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ Domain) +
  scale_fill_manual(values = c("Positive" = "#db4c34",
                               "Negative" = "#328de2",
                               "Not sig"  = "grey80")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Plant species", y = "Number of ASVs", fill = "Direction") +
  theme_de_bar

# ---- 5.2 Shared and unique ASVs per Domain (pseudo-venn) ----
asv_flag <- all_sig |>
  dplyr::distinct(Domain, ASV, Species_group) |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(names_from  = Species_group,
                     values_from = present,
                     values_fill = 0L) |>
  dplyr::mutate(Festuca = Festuca > 0, Lolium  = Lolium  > 0)

venn_counts <- asv_flag |>
  dplyr::group_by(Domain) |>
  dplyr::summarise(Festuca_only = sum(Festuca & !Lolium),
                   Lolium_only  = sum(!Festuca & Lolium),
                   Shared       = sum(Festuca & Lolium),
                   .groups      = "drop")

label_df <- venn_counts |>
  tidyr::pivot_longer(cols      = c(Festuca_only, Lolium_only, Shared),
                      names_to  = "Region", values_to = "Count") |>
  dplyr::mutate(x = dplyr::case_when(Region == "Festuca_only" ~ -1,
                                     Region == "Lolium_only"  ~  1,
                                     Region == "Shared"       ~  0),
                y = 0)

circle_df <- expand.grid(Domain = unique(asv_flag$Domain),
                         Species = c("Festuca", "Lolium"),
                         KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE) |>
  dplyr::mutate(x0 = ifelse(Species == "Festuca", -1, 1), y0 = 0, r = 1.5)

species_cols <- species_cols  # from utils_theme.R (Festuca / Lolium palette)

p_venn <- ggplot() +
  ggforce::geom_circle(data = circle_df,
                       aes(x0 = x0, y0 = y0, r = r, fill = Species),
                       alpha = 0.9, color = "black", linewidth = 0.6) +
  geom_text(data = label_df, aes(x = x, y = y, label = Count),
            fontface = "bold", size = 4.5, color = "white") +
  scale_fill_manual(values = species_cols, name = "Plant species") +
  coord_fixed(xlim = c(-3, 3), ylim = c(-2, 2), clip = "off") +
  facet_wrap(~ Domain) +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 11, margin = margin(t = 0.5, b = 1.5)),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold.italic"))

# -------------------------------------------------------------------
# 6) Shared ASVs: scatter of Festuca vs Lolium log2FC + top taxa table
# -------------------------------------------------------------------
# find ASVs present in both species within a domain
shared_ids <- all_sig |>
  dplyr::distinct(Domain, ASV, Species_group) |>
  dplyr::group_by(Domain, ASV) |>
  dplyr::summarise(n_species = dplyr::n_distinct(Species_group), .groups   = "drop") |>
  dplyr::filter(n_species == 2)

shared_df <- all_sig |>
  dplyr::inner_join(shared_ids, by = c("Domain", "ASV"))

shared_wide <- shared_df |>
  dplyr::select(ASV, Domain, Species_group, log2FoldChange,
                phylum, class, family, TrophicMode) |>
  tidyr::pivot_wider(names_from  = Species_group, values_from = log2FoldChange) |>
  dplyr::arrange(ASV) |>
  dplyr::arrange(Domain, dplyr::desc(Festuca)) |>
  dplyr::group_by(Domain) |>
  dplyr::mutate(ASV_short = paste0("ASV", dplyr::row_number())) |>
  dplyr::ungroup()

p_shared <- ggplot(shared_wide, aes(x = Festuca, y = Lolium, color = Domain)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_point(size = 4) +
  scale_color_manual(values = c("Bacteria" = "#1B9E77", "Fungi" = "#A6761D")) +
  labs(x = "log2 moisture-response (Festuca)", y = "log2 moisture-response (Lolium)",
       color = "Domain") +
  theme_classic() +
  theme(aspect.ratio    = 1,
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(color = species_cols["Festuca"], face = "bold", size = 12),
        axis.title.y = element_text(color = species_cols["Lolium"], face = "bold", size = 12))

# label points with ASV_short
p_shared_labeled <- p_shared +
  geom_text(aes(label = ASV_short), nudge_x = 0, nudge_y = 0.15,
            size = 3, fontface = "bold", check_overlap = TRUE)

# pick top 3 by Festuca log2FC (across domains)
top3 <- shared_wide |>
  dplyr::slice_max(Festuca, n = 3)

tax_info <- shared_df |>
  dplyr::distinct(ASV, Domain, phylum, class, family, TrophicMode)

summary_tab <- top3 |>
  dplyr::left_join(tax_info, by = c("ASV", "Domain")) |>
  dplyr::transmute(ASV = ASV_short, Domain, Phylum = phylum,
                   Class = class, Family = family, `Trophic mode` = TrophicMode)

# transpose to ASVs as columns for the table
asv_cols <- summary_tab$ASV
summary_tab_noASV <- summary_tab |>
  dplyr::select(-ASV)

summary_tab_T <- t(as.matrix(summary_tab_noASV)) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  tibble::rownames_to_column("Taxon")

colnames(summary_tab_T) <- c("Taxon", asv_cols)

theme_lines <- gridExtra::ttheme_minimal(core = list(
  fg_params = list(fontface = "bold", cex = 0.7),
  bg_params = list(col = "black", lwd = 0.5, fill = "white")),
  colhead = list(
    fg_params = list(fontface = "bold", cex = 0.8),
    bg_params = list(col = "black", lwd = 0.8, fill = "white")),
  rowhead = list(
    fg_params = list(fontface = "bold", cex = 0.7),
    bg_params = list(col = "black", lwd = 0.8, fill = "white")))

tab_grob <- gridExtra::tableGrob(summary_tab_T, rows = NULL, theme = theme_lines)

# combined panel: venn + scatter + table
p_shared_panel <- ggpubr::ggarrange(p_venn, NULL, p_shared_labeled, tab_grob,
                                    ncol = 1, heights = c(1, 0.15, 2, 1),
                                    labels = c("(a)", "", "(b)", "(c)"))

ggsave(file.path(results_dir, "Indicator_species_Moisture_overlap.png"), ## Figure S15
       p_shared_panel, width = 6.5, height = 8, units = "in", dpi = 600, bg = "white")

# -------------------------------------------------------------------
# 7) Lollipop plots by phylum (Bacteria & Fungi)
# -------------------------------------------------------------------

# split all_sig by domain
bac_sig_all <- dplyr::filter(all_sig, Domain == "Bacteria")
fun_sig_all <- dplyr::filter(all_sig, Domain == "Fungi")

theme_lolli <- theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "bottom",
        strip.background = element_rect(color = NA),
        strip.text = element_text(size = 11, color = "black", face = "bold.italic"),
        text = element_text(size = 11, face = "bold", color = "black"))

## --- Bacteria lollipop ---
phylum_counts_dir_bac <- bac_sig_all |>
  dplyr::count(Species_group, phylum, Diff.abund) |>
  tidyr::pivot_wider(names_from  = Diff.abund, values_from = n, values_fill = 0) |>
  dplyr::mutate(n_phylum = Positive + Negative)

bac_lolli <- bac_sig_all |>
  dplyr::left_join(phylum_counts_dir_bac, by = c("Species_group", "phylum")) |>
  dplyr::mutate(Phylum_label = paste0(
    phylum, "\n", "(Pos = ", Positive, ", Neg = ", Negative, " ASVs)")) |>
  dplyr::group_by(Species_group) |>
  dplyr::mutate(Phylum_label = forcats::fct_reorder(Phylum_label, n_phylum)) |>
  dplyr::ungroup()

p_lolli_bac <- ggplot(bac_lolli, aes(x = log2FoldChange, y = Phylum_label, color = phylum)) +
  geom_segment(aes(x = 0, xend = log2FoldChange, y = Phylum_label, yend = Phylum_label),
               linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 2, position = position_jitter(height = 0.15)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ Species_group, scales = "free_y") +
  labs(x = "log2 fold-change along moisture gradient",
       y = "Phylum (number of moisture-responsive ASVs)", color = "Bacterial phylum") +
  theme_lolli

# add legend labels with total ASVs per bacterial phylum
phylum_totals_bac <- tax_Bac |>
  dplyr::count(phylum, name = "Total_ASVs_phylum")

legend_labels_bac <- phylum_totals_bac |>
  dplyr::mutate(phylum_label = paste0(phylum, " (", Total_ASVs_phylum, ")")) |>
  dplyr::select(phylum, phylum_label)

p_lolli_bac <- p_lolli_bac +
  scale_color_discrete(labels = setNames(legend_labels_bac$phylum_label,
                                         legend_labels_bac$phylum))

## --- Fungi lollipop ---
phylum_counts_dir_fun <- fun_sig_all |>
  dplyr::count(Species_group, phylum, Diff.abund) |>
  tidyr::pivot_wider(names_from  = Diff.abund, values_from = n, values_fill = 0) |>
  dplyr::mutate(n_phylum = Positive + Negative)

fun_lolli <- fun_sig_all |>
  dplyr::left_join(phylum_counts_dir_fun, by = c("Species_group", "phylum")) |>
  dplyr::mutate(Phylum_label = paste0(
      phylum, "\n", "(Pos = ", Positive, ", Neg = ", Negative, " ASVs)")) |>
  dplyr::group_by(Species_group) |>
  dplyr::mutate(Phylum_label = forcats::fct_reorder(Phylum_label, n_phylum)) |>
  dplyr::ungroup() |>
  dplyr::arrange(ASV) |>
  dplyr::mutate(ASV_short = paste0("ASV", dplyr::row_number()))

p_lolli_fun <- ggplot(fun_lolli, aes(x = log2FoldChange, y = Phylum_label, color = phylum)) +
  geom_segment(aes(x = 0, xend = log2FoldChange, y = Phylum_label, yend = Phylum_label),
               linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 2, position = position_jitter(height = 0.15)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ Species_group, scales = "free_y") +
  labs(x = "log2 fold-change along moisture gradient",
       y = "Phylum (number of moisture-responsive ASVs)", color = "Fungal phylum") +
  theme_lolli

phylum_totals_fun <- tax_Fun |>
  dplyr::count(phylum, name = "Total_ASVs_phylum")

legend_labels_fun <- phylum_totals_fun |>
  dplyr::mutate(phylum_label = paste0(phylum, " (", Total_ASVs_phylum, ")")) |>
  dplyr::select(phylum, phylum_label)

p_lolli_fun <- p_lolli_fun +
  scale_color_discrete(labels = setNames(legend_labels_fun$phylum_label,
                                         legend_labels_fun$phylum))

# annotate a few fungal ASVs with text boxes
annot_ids <- c("ASV46", "ASV33", "ASV4")
annot_df  <- fun_lolli |>
  dplyr::filter(ASV_short %in% annot_ids)

p_lolli_fun_labeled <- p_lolli_fun +
  geom_text(data = annot_df, aes(label = ASV_short), nudge_x = -0.15, nudge_y = 0.05,
            size = 3, fontface = "bold", color = "black")

# three small info boxes (labels to ASVs)
p_textbox1 <- ggplot() +
  ggtext::geom_richtext(aes(x = 0.5, y = 0.5, label = "<b>ASV46</b><br
  Class: <b>Tremellomycetes</b><br
  Family: <b>Trichosporonaceae</b><br>
  Trophic mode: <b>Saprotroph</b>"),
                        fill = alpha("gray95", 0.9), size = 3,
                        hjust = 0.5, vjust = 0.5, lineheight = 1.5) +
  theme_void()

p_textbox2 <- ggplot() +
  ggtext::geom_richtext(aes(x = 0.5, y = 0.5, label = "<b>ASV33</b><br>
  Class: <b>Dothideomycetes</b><br>
  Trophic mode: <b>Saprotroph–Symbiotroph</b>"),
                        fill = alpha("gray95", 0.9), size = 3,
                        hjust = 0.5, vjust = 0.5, lineheight = 1.5) +
  theme_void()

p_textbox3 <- ggplot() +
  ggtext::geom_richtext(aes(x = 0.5, y = 0.5, label = "<b>ASV4</b><br>
  Class: <b>Sordariomycetes</b><br>
  Family: <b>Hypocreaceae</b><br>
  Trophic mode: <b>Pathotroph–Saprotroph</b>"),
                        fill = alpha("gray95", 0.9), size = 3,
                        hjust = 0.5, vjust = 0.5, lineheight = 1.5) +
  theme_void()

p_lolli_fun_annot <- cowplot::ggdraw(p_lolli_fun_labeled) +
  cowplot::draw_plot(p_textbox1, x = 0.24, y = 0.35, width = 0.35, height = 0.25) +
  cowplot::draw_plot(p_textbox2, x = 0.16, y = 0.56, width = 0.45, height = 0.25) +
  cowplot::draw_plot(p_textbox3, x = 0.51, y = 0.48, width = 0.45, height = 0.25)

# combined lollipop figure (Bacteria + Fungi)
pc1 <- ggpubr::ggarrange(p_lolli_bac, p_lolli_fun_annot, ncol = 1,
                         heights = c(1.5, 1), labels = c("(a)", "(b)"),
                         font.label = list(size = 12))

ggsave(file.path(results_dir, "Indicator_species_Moisture_taxa.png"), ## Figure 4
       pc1, width = 14, height = 12, units = "in", dpi = 600, bg = "white")


