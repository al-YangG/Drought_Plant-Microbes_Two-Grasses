##-## Drought Responses in Two Grassland Plant-Microbe Systems ##-##
##-## Author: Gang Yang ##-##
##-## Last edit: December 2025 ##-##
##
## R/utils_indicator_taxa.R

# Example helper to count significant indicators per group
count_signif_indicators <- function(ind_table,
                                    p_col = "p_value",
                                    alpha = 0.05,
                                    group_cols = c("Domain", "Compartment", "Group")) {
  ind_table |>
    dplyr::filter(.data[[p_col]] <= alpha) |>
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "n_indicators")
}

