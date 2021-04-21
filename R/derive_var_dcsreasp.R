
derive_var_dcsreasp <- function(dataset,
                              dataset_ds = ds,
                              filter_ds = exprs(DSCAT == "DISPOSITION EVENT" & DSSCAT =="STUDY COMPLETION/EARLY DISCONTINUATION")
){

  derive_merged_vars(dataset,
                     dataset_add = dataset_ds,
                     filter_add = filter_ds,
                     by_vars = exprs(USUBJID),
                     new_vars = exprs(temp___= DSTERM)
  )%>%
    mutate(DCSREASP = case_when(
      temp___ != "COMPLETED" ~ temp___,
      TRUE ~ ""
    )
    ) %>%
    select(-ends_with("___"))

}

