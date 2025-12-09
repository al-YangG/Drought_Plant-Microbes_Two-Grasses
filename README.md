# Drought Responses in Two Grassland Plant-Microbe Systems Under Climate Change

This repository contains the full analysis workflow and data used in the study:

**"Contrasting drought responses in two grassland plant–microbe systems under climate change."**

The project examines how two grass species with contrasting ecological strategies — **Festuca rubra** (stress-tolerant) and **Lolium perenne** (competitive) — respond to a continuous moisture gradient under ambient vs. future climate conditions (warming + elevated CO2). Analyses integrate plant performance, traits, microbial diversity, composition, functional attributes, and structural equation modelling.

---

## Repository Structure

```
Drought_Plant-Microbes_Two-Grasses/
│
├── R/    # All R scripts (analysis workflow)
│     0_setup_paths.R
│     1_plant_performance.R
│     1b_plant_performance_visualization.R
│     2_phyloseq_objects.R
│     2b_qc_venn_rarefaction.R
│     3_alpha_diversity_models.R
│     3b_alpha_diversity_visualization.R
│     4_ordination_RDA_PCA.R
│     5_variance_explained.R
│     6_taxonomic_composition.R
│     6b_taxonomic_top5_models.R
│     7_indicator_species_indicspecies.R
│     7b_indicator_species_DESeq2_moisture.R
│     7c_indicator_DESeq2_summaries_vs_indicspecies.R
│     8_functional_pathways_Tax4Fun2.R
│     8b_fungal_functional_groups.R
│     9_structural_equation_models.R
│     utils_packages.R
│     utils_theme.R
│     utils_phyloseq.R
│     utils_stat.R
│     utils_indicator_taxa.R
│
├── Data/
│     1_Plant performance/
│     2_Phyloseq objects/
│     3_Alpha diversity/
│     4_PCA scores/
│     5_Variation explained/
│     6_Taxonomic composition/
│     7_Indicator species/
│     8_Functional group/
│     9_SEM/
│
├── results/    # Output figures and tables (created when scripts run)
│
├── .gitignore
├── Drought_Plant-Microbes_Two-Grasses.Rproj
└── README.md
```

---

## How to Run the Project

### 1. Open the project
Double-click:
```
Drought_Plant-Microbes_Two-Grasses.Rproj
```

### 2. Install and load required packages
All packages are handled by:
```r
source("R/utils_packages.R")
load_project_packages()
```

### 3. Run scripts in order
Scripts are numbered following the analysis workflow:

- 0_setup_paths.R — project setup  
- 1 / 1b — plant performance  
- 2 / 2b — phyloseq objects, QC, rarefaction  
- 3 / 3b — alpha diversity analyses  
- 4 — ordination (PCA, RDA)  
- 5 — variation partitioning  
- 6 / 6b — taxonomic composition and top taxa  
- 7 / 7b / 7c — indicator taxa (indicspecies + DESeq2)  
- 8 / 8b — predicted functional pathways  
- 9 — structural equation models  

Each script writes outputs to the `results/` folder.

---

## Analyses Included

- Plant performance and trait responses  
- Microbial richness and diversity  
- Community composition (bacteria and fungi)  
- Indicator species (indicspecies and DESeq2)  
- Functional prediction (Tax4Fun2, fungal guilds)  
- Variation explained  
- Species-specific structural equation models  
- Correlation networks and integrated pathway analysis  

---

## Data Availability

The `Data/` directory contains analysis-ready datasets.  
Large raw sequencing files are not stored in this repository and will be provided separately upon publication.

---

## Reproducibility Notes

- Paths use `here::here()` or relative paths  
- Scripts are fully automated  
- Output is written to `results/`  
- R version ≥ 4.2 recommended  

---

## Contact

For questions or suggestions, please open a GitHub Issue or contact the corresponding author: yangg@natur.cuni.cz.
