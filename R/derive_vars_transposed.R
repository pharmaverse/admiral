derive_vars_transposed <- function(dataset,
                                   dataset_merge,
                                   by_vars,
                                   key_var,
                                   value_var,
                                   filter = NULL) {
  key_var <- assert_symbol(enquo(key_var))
  value_var <- assert_symbol(enquo(value_var))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(dataset_merge, required_vars = quo_c(by_vars, key_var, value_var))

  dataset_transposed <- dataset_merge %>%
    filter_if(filter) %>%
    spread(key = !!key_var, value = !!value_var)

  left_join(dataset, dataset_transposed, by = vars2chr(by_vars))
}

derive_vars_atc <- function(dataset,
                            dataset_facm,
                            by_vars = vars(USUBJID, CMREFID = FAREFID)) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(dataset_facm, required_vars = vars(!!!by_vars, FAGRPID, FATESTCD, FASTRESC))

  dataset %>%
    derive_vars_transposed(
      select(dataset_facm, !!!unname(by_vars), FAGRPID, FATESTCD, FASTRESC),
      by_vars = by_vars,
      key_var = FATESTCD,
      value_var = FASTRESC,
      filter = str_detect(FATESTCD, "^CMATC[1-4](CD)?$")
    ) %>%
    select(-starts_with("FA")) %>%
    rename_at(vars(starts_with("CMATC")), ~str_remove(.x, "^CM"))
}
