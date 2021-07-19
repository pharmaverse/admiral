derive_derived_param <- function(dataset,
                                 filter = NULL,
                                 parameters,
                                 by_vars,
                                 analysis_value,
                                 set_values_to
) {
  # checking and quoting
  assert_vars(by_vars)
  assert_data_frame(dataset,
                    required_vars = vars(!!!by_vars, AVAL))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_character_vector(parameters)

  # select observations and variables required for new observations
  data <- dataset %>%
    filter_if(filter) %>%
    filter(PARAMCD %in% parameters)

  keep_vars <- get_constant_vars(data, by_vars = by_vars)
  data <- data %>%
    select(c(vars2chr(keep_vars),
             "PARAMCD",
             "AVAL"))

  signal_duplicate_records(
    data,
    by_vars = vars(!!!by_vars, PARAMCD),
    msg = paste("The filtered input dataset contains duplicate records with respect to",
                enumerate(c(vars2chr(by_vars), "PARAMCD")),
                "\nPlease ensure that the variables specified for `by_vars` and `PARAMCD`",
                "are a unique key of the input data set restricted by the condition",
                "specified for `filter` and to the parameters specified for `parameters`.")
  )

  # horizontalize data, AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  hori_data <- data %>%
    spread(key = PARAMCD,
           value = AVAL,
           sep = ".")
  names(hori_data) <- map_chr(names(hori_data), str_replace, "PARAMCD.", "AVAL.")

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  hori_data <- hori_data %>%
    mutate(AVAL = !!enquo(analysis_value),
           !!!set_values_to) %>%
    select(-starts_with("AVAL."))

  bind_rows(dataset, hori_data)
}
