##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## R/utils_phyloseq.R

# Build a phyloseq object from three CSVs
build_phyloseq_from_csv <- function(asv_fp, taxa_fp, meta_fp,
                                    taxa_are_rows = TRUE) {
  asv  <- read.csv(asv_fp, row.names = 1, check.names = FALSE)
  taxa <- read.csv(taxa_fp, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_fp, row.names = 1, check.names = FALSE)
  
  phyloseq::phyloseq(phyloseq::sample_data(meta),
                     phyloseq::otu_table(as.matrix(asv), taxa_are_rows = taxa_are_rows),
                     phyloseq::tax_table(as.matrix(taxa)))
}

# Rarefy to min or to median if zeros present
rarefy_ps_safe <- function(ps, seed = 100) {
  set.seed(seed)
  ss_min <- min(phyloseq::sample_sums(ps))
  if (ss_min == 0) {
    message("Some samples have 0 reads, using median depth instead.")
    target <- stats::median(phyloseq::sample_sums(ps))
  } else {
    target <- ss_min
  }
  ps_r <- phyloseq::rarefy_even_depth(ps, sample.size = target)
  ps_r <- phyloseq::filter_taxa(ps_r, function(x) sum(x) > 0, TRUE)
  ps_r
}

# Alpha-diversity table from a rarefied phyloseq object
compute_alpha_div <- function(ps) {
  count <- as.data.frame(phyloseq::otu_table(ps))
  x     <- t(count)
  
  Shannon       <- vegan::diversity(x)
  Inv_Simpson   <- vegan::diversity(x, index = "invsimpson")
  S             <- vegan::specnumber(x)
  Pielou_even   <- Shannon / log(S)
  Simpson_even  <- Inv_Simpson / S
  est           <- vegan::estimateR(x)
  Richness      <- est[1, ]
  Chao1         <- est[2, ]
  ACE           <- est[4, ]
  
  cbind(Shannon           = Shannon,
        Inv_Simpson       = Inv_Simpson,
        Pielou_evenness   = Pielou_even,
        Simpson_evenness  = Simpson_even,
        Richness          = Richness,
        Chao1             = Chao1,
        ACE               = ACE)
}
