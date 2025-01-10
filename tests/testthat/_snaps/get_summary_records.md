# get_summary_records Test 1: Message sent to users

    Code
      df <- input %>% get_summary_records(by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(AVAL = mean(AVAL, na.rm = TRUE), DTYPE = "AVERAGE")) %>%
        dplyr::mutate(AVAL = round(AVAL))
    Message
      `get_summary_records()` was deprecated in admiral 1.2.0.
      i Please use `derive_summary_records()` instead.
      x This message will turn into a warning for at least one year.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation

