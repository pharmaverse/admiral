
derive_var_dcsreas <- function(dataset,
                               dataset_ds,
                               filter_ds = exprs(DSCAT == "DISPOSITION EVENT" & DSSCAT == "STUDY COMPLETION/EARLY DISCONTINUATION")) {
  derive_merged_vars(
    dataset,
    dataset_add = dataset_ds,
    filter_add = filter_ds,
    new_vars = exprs(temp___ = DSDECOD)
  ) %>%
    mutate(DCSREAS = case_when(
      temp___ != "COMPLETED" ~ temp___,
      TRUE ~ NA_character_
    )) %>%
    select(-ends_with("___"))
}
