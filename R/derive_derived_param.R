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
  if (!quo_is_null(filter)) {
    data <- dataset %>% filter(!!filter)
  }
  else {
    data <- dataset
  }
  data <- data %>%
    select(!!!by_vars, PARAMCD, AVAL) %>%
    filter(PARAMCD %in% parameters)

  # horizontalize data, AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  hori_data <- data %>%
    spread(key = PARAMCD,
           value = AVAL,
           sep = ".")
  names(hori_data) <- map_chr(names(hori_data), str_replace, "PARAMCD.", "AVAL.")

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  hori_data %>%
    mutate(AVAL = !!enquo(analysis_value),
           !!!set_values_to) %>%
    select(-starts_with("AVAL."))
}
