# derive_param_tte Test 6: an error is issued if some of the by variables are missing

    Code
      derive_param_tte(dataset_adsl = adsl, by_vars = exprs(AEBODSYS, AEDECOD),
      start_date = TRTSDT, event_conditions = list(ttae), censor_conditions = list(
        eos), source_datasets = list(adsl = adsl, ae = ae), set_values_to = exprs(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))), PARAM = paste(
          "Time to First", AEDECOD, "Adverse Event"), PARCAT1 = "TTAE", PARCAT2 = AEDECOD))
    Condition
      Error in `signal_duplicate_records()`:
      ! Required variables `AEBODSYS` and `AEDECOD` are missing in `dataset`

# derive_param_tte Test 7: errors if all by vars are missing in all source datasets

    Code
      derive_param_tte(dataset_adsl = adsl, by_vars = exprs(AEBODSYS), start_date = TRTSDT,
      event_conditions = list(ttae), censor_conditions = list(eos), source_datasets = list(
        adsl = adsl, ae = ae), set_values_to = exprs(PARAMCD = paste0("TTAE",
        as.numeric(as.factor(AEDECOD))), PARAM = paste("Time to First", AEDECOD,
        "Adverse Event"), PARCAT1 = "TTAE", PARCAT2 = AEDECOD))
    Condition
      Error in `signal_duplicate_records()`:
      ! Required variable `AEBODSYS` is missing in `dataset`

# derive_param_tte Test 8: errors if PARAMCD and by_vars are not one to one

    Code
      derive_param_tte(dataset_adsl = adsl, by_vars = exprs(AEDECOD), start_date = TRTSDT,
      event_conditions = list(ttae), censor_conditions = list(eos), source_datasets = list(
        adsl = adsl, ae = ae), set_values_to = exprs(PARAMCD = "TTAE", PARCAT2 = AEDECOD))
    Condition
      Error in `signal_duplicate_records()`:
      ! Required variable `AEDECOD` is missing in `dataset`

# derive_param_tte Test 9: errors if set_values_to contains invalid expressions

    Code
      derive_param_tte(dataset_adsl = adsl, by_vars = exprs(AEDECOD), start_date = TRTSDT,
      event_conditions = list(ttae), censor_conditions = list(eos), source_datasets = list(
        adsl = adsl, ae = ae), set_values_to = exprs(PARAMCD = paste0("TTAE",
        as.numeric(as.factor(AEDECOD))), PARAM = past("Time to First", AEDECOD,
        "Adverse Event"), PARCAT1 = "TTAE", PARCAT2 = AEDECOD))
    Condition
      Error in `signal_duplicate_records()`:
      ! Required variable `AEDECOD` is missing in `dataset`

# list_tte_source_objects Test 13: error is issued if package does not exist

    Code
      list_tte_source_objects(package = "tte")
    Condition
      Error in `list_tte_source_objects()`:
      ! No package called tte is installed and hence no <tte_source> objects are available.

