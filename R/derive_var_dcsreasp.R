
derive_var_dcsreasp <- function(dataset,
                                dataset_ds,
                                filter_ds = exprs(DSCAT == "DISPOSITION EVENT" & DSSCAT == "STUDY COMPLETION/EARLY DISCONTINUATION")) {
  derive_merged_vars(
    dataset,
    dataset_add = dataset_ds,
    filter_add = filter_ds,
    new_vars = exprs(term___ = DSTERM, decod___ = DSDECOD)
  ) %>%
    mutate(
      DCSREASP = case_when(
        decod___ != "COMPLETED" ~ term___,
        TRUE ~ NA_character_
      )
    ) %>%
    select(-ends_with("___"))
}
