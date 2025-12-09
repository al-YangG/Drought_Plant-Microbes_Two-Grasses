##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## scripts/8b_fungal_functional_groups.R
##
## Tasks:
##  1) Add fungal functional groups (TrophicMode / Guild) to ps-Fun.
##  2) Summarise trophic modes across Compartment, Moisture level (M1–M6),
##     Species, and Climate, and make stacked barplots.
##  3) Fit beta-regression models to test effects of experimental factors
##     on trophic-mode relative abundances (global model).
##  4) Focus on Compartment effects: Mesocosm-level means, pairwise
##     comparisons, and barplots with significance letters.
## -------------------------------------------------------------------

# -------------------------------------------------------------------
# 0) Setup
# -------------------------------------------------------------------
library(phyloseq)
library(dplyr)
library(forcats)
library(scales)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(glmmTMB)
library(car)
library(emmeans)
library(multcomp)
library(tidyr)
library(purrr)

root        <- here::here()
ps_dir      <- file.path(root, "Data", "2_Phyloseq objects", "Output", "ps_new")
fg_dir      <- file.path(root, "Data", "8_Functional group")
results_dir <- file.path(root, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# 1) Add functional group to fungal phyloseq object
# -------------------------------------------------------------------
ps_fun <- readRDS(file.path(ps_dir, "ps-Fun.rds"))

# Helper: add MoistureLevel (M1–M6) to sample_data
add_M_levels <- function(ps) {
  sd <- as.data.frame(sample_data(ps))
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
  sample_data(ps)$MoistureLevel <- sd$MoistureLevel
  ps
}

ps_fun <- add_M_levels(ps_fun)

# Read functional table and align to ASVs
fg_raw <- read.csv(file.path(fg_dir, "fungal guild table.csv"))
fg <- fg_raw %>%
  select(ASVid, Group, TrophicMode, Guild)

fg2 <- fg %>%
  filter(ASVid %in% taxa_names(ps_fun)) %>%
  distinct(ASVid, .keep_all = TRUE) %>%
  slice(match(taxa_names(ps_fun), ASVid))

tt <- as.data.frame(tax_table(ps_fun), stringsAsFactors = FALSE)

for (col in c("Group", "TrophicMode", "Guild")) {
  if (!col %in% colnames(tt)) tt[[col]] <- NA_character_
}

matched_idx <- which(!is.na(fg2$ASVid))
tt$Group[matched_idx]       <- fg2$Group[matched_idx]
tt$TrophicMode[matched_idx] <- fg2$TrophicMode[matched_idx]
tt$Guild[matched_idx]       <- fg2$Guild[matched_idx]

tax_table(ps_fun) <- as.matrix(tt)

message("Matched functional annotations for ",
        sum(!is.na(fg2$ASVid)), " of ", ntaxa(ps_fun), " ASVs.")

# Save updated object
saveRDS(ps_fun, file.path(ps_dir, "ps-Fun_withGuilds.rds"))

# -------------------------------------------------------------------
# 2) Relative abundance of trophic modes across factors
# -------------------------------------------------------------------
ps_fun <- readRDS(file.path(ps_dir, "ps-Fun_withGuilds.rds"))

# Explicit "Unassigned"
tt <- as.data.frame(tax_table(ps_fun), stringsAsFactors = FALSE)
tt$TrophicMode[is.na(tt$TrophicMode) | tt$TrophicMode == ""] <- "Unassigned"
tax_table(ps_fun) <- as.matrix(tt)

# Aggregate by TrophicMode
ps_fun_tm    <- tax_glom(ps_fun, taxrank = "TrophicMode", NArm = FALSE)
ps_fun_tm_rel <- transform_sample_counts(
  ps_fun_tm,
  function(x) if (sum(x) > 0) x / sum(x) else x
)

df_tm <- psmelt(ps_fun_tm_rel)

# Desired TrophicMode order
mode_order <- c("Pathotroph", "Saprotroph", "Symbiotroph", "Pathotroph-Saprotroph",
                "Pathotroph-Symbiotroph", "Saprotroph-Symbiotroph",
                "Pathotroph-Saprotroph-Symbiotroph", "Unassigned")

# --- Summaries per factor ---

# Compartment (rename to endosphere / rhizosphere)
tm_comp <- df_tm %>%
  group_by(Compartment, TrophicMode) %>%
  summarise(mean_RA = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Compartment = dplyr::recode(Compartment,
                                     "Leaf" = "Leaf endosphere",
                                     "Root" = "Root endosphere",
                                     "Soil" = "Rhizosphere"),
    Compartment = factor(Compartment,
                         levels = c("Rhizosphere", "Root endosphere", "Leaf endosphere"))) %>%
  group_by(Compartment) %>%
  mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TrophicMode = factor(TrophicMode, levels = mode_order))

# Moisture levels (M1–M6)
tm_moist <- df_tm %>%
  group_by(MoistureLevel, TrophicMode) %>%
  summarise(mean_RA = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(MoistureLevel = factor(MoistureLevel, levels = paste0("M", 6:1))) %>%
  group_by(MoistureLevel) %>%
  mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TrophicMode = factor(TrophicMode, levels = mode_order))

# Species
tm_species <- df_tm %>%
  filter(!is.na(Species)) %>%
  group_by(Species, TrophicMode) %>%
  summarise(mean_RA = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Species = factor(Species, levels = c("Lolium", "Festuca"))) %>%
  group_by(Species) %>%
  mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TrophicMode = factor(TrophicMode, levels = mode_order))

# Climate
tm_climate <- df_tm %>%
  group_by(Climate, TrophicMode) %>%
  summarise(mean_RA = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Climate = factor(Climate, levels = c("Current", "Future"))) %>%
  group_by(Climate) %>%
  mutate(mean_RA = mean_RA / sum(mean_RA, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TrophicMode = factor(TrophicMode, levels = mode_order))

# Combine for a single faceted plot
comb_tm <- bind_rows(
  tm_comp    %>% transmute(GroupType = "Compartment",    X = Compartment,    TrophicMode, mean_RA),
  tm_moist   %>% transmute(GroupType = "Moisture level", X = MoistureLevel, TrophicMode, mean_RA),
  tm_species %>% transmute(GroupType = "Species",        X = Species,       TrophicMode, mean_RA),
  tm_climate %>% transmute(GroupType = "Climate",        X = Climate,       TrophicMode, mean_RA)
) %>%
  mutate(GroupType = factor(GroupType,
                            levels = c("Compartment", "Moisture level", "Species", "Climate"))) %>%
  group_by(GroupType) %>%
  mutate(X          = forcats::fct_rev(factor(X, levels = unique(X))),
         TrophicMode = forcats::fct_rev(TrophicMode)) %>%
  ungroup()

tm_levels <- levels(forcats::fct_rev(factor(comb_tm$TrophicMode)))
pal_fun_tm <- colorRampPalette(brewer.pal(8, "Set2"))(length(tm_levels))

theme_tm_overview <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 12, color = "black"),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        panel.spacing = unit(0.5, "lines"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

p_fun_tm_vert <- ggplot(comb_tm, aes(y = X, x = mean_RA, fill = TrophicMode)) +
  geom_col(width = 0.7, position = "stack", color = NA) +
  scale_fill_manual(values = setNames(rev(pal_fun_tm), rev(tm_levels)), drop = FALSE,
                    guide = guide_legend(reverse = TRUE, nrow = 3, byrow = TRUE)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Mean relative abundance (%)", y = NULL, fill = "Trophic mode",
       title = "Fungal trophic modes across environmental and host factors") +
  facet_grid(rows = vars(GroupType), scales = "free_y", space = "free_y") +
  theme_tm_overview

ggsave(file.path(results_dir, "Fungal_trophic_modes_overview.png"),
       p_fun_tm_vert, width = 10, height = 12, units = "in", bg = "white", dpi = 600)

# -------------------------------------------------------------------
# 3) Beta regression: global model for each trophic mode
# -------------------------------------------------------------------
df_tm_beta <- psmelt(ps_fun_tm_rel) %>%
  select(Sample, Abundance, TrophicMode, Compartment,
         Moisture, Climate, Species, Mesocosm) %>%
  mutate(Moisture_sc = as.numeric(scale(Moisture)),
         Mesocosm    = factor(Mesocosm))

# helper: beta regression by taxon with generic formula
run_beta_by_taxon <- function(df, tax_col, fixed_formula) {
  split_df <- split(droplevels(df), df[[tax_col]])
  
  purrr::map_df(split_df, function(sub) {
    n <- nrow(sub)
    sub$Abund01 <- (sub$Abundance * (n - 1) + 0.5) / n
    
    ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
    
    m <- glmmTMB(formula = stats::as.formula(fixed_formula),
                 data    = sub,
                 family  = beta_family(link = "logit"),
                 control = ctrl)
    
    a <- car::Anova(m, type = 3)
    
    data.frame(Taxon = unique(sub[[tax_col]]),
               Effect = rownames(a),
               Stat   = a$Chisq,
               P      = a$`Pr(>Chisq)`,
               Test   = "Wald χ² (beta)",
               row.names = NULL)
  })
}

options(contrasts = c("contr.sum", "contr.poly"))

res_tm <- run_beta_by_taxon(df    = df_tm_beta, tax_col = "TrophicMode",
                            fixed_formula = "Abund01 ~ Compartment * Moisture_sc * Climate * Species + (1 | Mesocosm)")

write.csv(res_tm, row.names = FALSE,
          file.path(results_dir, "TrophicMode_effects_global.csv")) ## Table S11

# -------------------------------------------------------------------
# 4) Compartment-only model & barplots with letters
# -------------------------------------------------------------------
comp_order <- c("Leaf endosphere", "Root endosphere", "Rhizosphere")
tm_order   <- mode_order[mode_order != "Unassigned"]

df_tm_comp <- df_tm_beta %>%
  mutate(Compartment = dplyr::recode(Compartment,
                                     "Leaf" = "Leaf endosphere",
                                     "Root" = "Root endosphere",
                                     "Soil" = "Rhizosphere"),
         Compartment = factor(Compartment, levels = comp_order),
         TrophicMode = as.character(TrophicMode))

# Mesocosm means to avoid pseudoreplication
df_tm_ms <- df_tm_comp %>%
  group_by(Mesocosm, Compartment, TrophicMode) %>%
  summarise(Abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

# mean ± SE across mesocosms
sum_ms <- df_tm_ms %>%
  filter(TrophicMode != "Unassigned") %>%
  mutate(TrophicMode = factor(TrophicMode, levels = tm_order)) %>%
  group_by(Compartment, TrophicMode) %>%
  summarise(n = n(), mean = mean(Abund, na.rm = TRUE),
            se = sd(Abund,  na.rm = TRUE) / sqrt(n), .groups = "drop")

# significance letters per trophic mode
get_letters_one_mode <- function(sub_df) {
  n <- nrow(sub_df)
  sub_df$Abund01 <- (sub_df$Abundance * (n - 1) + 0.5) / n
  
  ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
  m <- glmmTMB(Abund01 ~ Compartment + (1 | Mesocosm), data = sub_df,
               family = beta_family(link = "logit"), control = ctrl)
  
  emm     <- emmeans(m, ~ Compartment, type = "response")
  cld_tbl <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
  
  as.data.frame(cld_tbl) %>%
    transmute(Compartment, .group = gsub("\\s+", "", .group))
}

letters_df <- df_tm_comp %>%
  filter(TrophicMode != "Unassigned") %>%
  mutate(TrophicMode = factor(TrophicMode, levels = tm_order)) %>%
  group_split(TrophicMode, keep = TRUE) %>%
  map_df(function(x) {
    tm <- as.character(unique(x$TrophicMode))
    res <- tryCatch(get_letters_one_mode(x), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    res$TrophicMode <- tm
    res
  }) %>%
  mutate(TrophicMode = factor(TrophicMode, levels = tm_order),
         Compartment = factor(Compartment, levels = comp_order))

sum_ms_letters <- sum_ms %>%
  left_join(letters_df, by = c("TrophicMode", "Compartment")) %>%
  mutate(label_y = (mean + se) * 100 + 0.03)

theme_bar <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        legend.position = c(0.9, 0.05),
        legend.justification = c("right", "bottom"),
        legend.key.size  = unit(0.6, "cm"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        panel.spacing = unit(0.5, "lines"))

p_bar <- ggplot(sum_ms, aes(Compartment, mean * 100, fill = Compartment)) +
  geom_col(width = 0.7, color = "grey20") +
  geom_errorbar(aes(ymin = pmax((mean - se) * 100, 0),
                    ymax = pmin((mean + se) * 100, 100)), width = 0.15) +
  facet_wrap(~ TrophicMode, scales = "free_y", ncol = 2) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.15))) +
  labs(y = "Relative abundance (%)", x = NULL) +
  theme_bar

p_bar_letters <- p_bar +
  geom_text(data = sum_ms_letters %>% filter(!is.na(.group)),
            aes(x = Compartment, y = label_y, label = .group),
            size = 3.5, fontface = "bold")

p_combined <- ggarrange(p_fun_tm_vert, p_bar_letters,
                        nrow = 1, ncol = 2, labels = c("(a)", "(b)"))

ggsave(file.path(results_dir, "Fungal_functional_groups.png"), ## Figure S13
       p_combined, width = 16, height = 10, units = "in", bg = "white", dpi = 600)

