
derive_var_eosstt <- function(dataset,
                              dataset_ds,
                              filter_ds = expr(DSCAT == "DISPOSITION EVENT" & DSSCAT == "STUDY COMPLETION/EARLY DISCONTINUATION")) {
  derive_merged_vars(
    dataset,
    dataset_add = dataset_ds,
    filter_add = filter_ds,
    new_vars = exprs(temp___ = DSDECOD)
  ) %>%
    mutate(EOSSTT = case_when(
      temp___ == "COMPLETED" ~ "COMPLETED",
      temp___ != "COMPLETED" ~ "DISCONTINUED",
      TRUE ~ "ONGOING"
    )) %>%
    select(-ends_with("___"))
}
