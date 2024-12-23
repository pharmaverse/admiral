# get_summary_records Test 1: Summarize average of triplicate ECG interval values

    Code
      df <- input %>% get_summary_records(by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(AVAL = mean(AVAL, na.rm = TRUE), DTYPE = "AVERAGE")) %>%
        dplyr::mutate(AVAL = round(AVAL))
    Message
      `get_summary_records()` was deprecated in admiral 1.2.0.
      i Please use `derive_summary_records()` instead.
      x This message will turn into a warning with release of 1.3.0
      i https://pharmaverse.github.io/admiral/reference/derive_summary_records.html

