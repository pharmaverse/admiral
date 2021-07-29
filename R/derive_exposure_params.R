derive_exposure_params <- function(dataset,
                                   by_vars,
                                   new_param ,
                                   input_param  ,
                                   fns,
                                   filter_rows = NULL,
                                   set_values_to = NULL,
                                   drop_values_from = NULL
) {

  assert_data_frame(dataset,
                    required_vars = quo_c(by_vars, vars( PARAMCD, AVAL, AVALC,  ASTDTM, AENDTM)))
  assert_character_scalar(new_param)
  assert_character_scalar(input_param)
  #HOW to do this?
  #assert_fns???
  assert_vars(drop_values_from, optional = TRUE)
  by_vars <- assert_vars(by_vars)
  filter_rows <- assert_filter_cond(enquo(filter_rows), optional = TRUE)

  if (!is.null(set_values_to)) {
    assert_that(is_quosures(set_values_to),
                msg = str_glue("`set_values_to` must be a `vars()` object, \\
                             not {friendly_type(typeof(set_values_to))}.")
    )
  }
  #TO DO
  #Assert drop_values from

  #TO DO
  #add this when available
  #assert_paramcd_does_not_exist(dataset, new_param)
  params_available <- unique(dataset$PARAMCD)
  assert_character_vector(input_param,
                          values = params_available)

  subset_ds <- dataset %>%
    filter_if(filter_rows)

  add_data <- subset_ds %>%
    filter(PARAMCD == input_param) %>%
    derive_summary_records(
      by_vars = by_vars,
      fns = fns,
      set_values_to = vars(PARAMCD___ = !!new_param),
      drop_values_from = drop_values_from
    ) %>%
    filter(PARAMCD___ == new_param) %>%
    mutate(PARAMCD = PARAMCD___) %>%
    #TO UPDATE
    # not sure why the derive_summary_records() fns render a AVAL.x and AVAL.y...
    # when i compute summary for TPDOSE
    #If i remove the summary for TPDOSE it works fine...
    # AVAL.x retain the value of the input param, AVAL.y has the correct result
    select(-ends_with("___"), -ends_with(".x")) #%>%
    mutate(
      AVAL = coalesce(!!!select(
        ., starts_with("AVAL"),
        -ends_with("C"),
        -ends_with("C.y")
      )),
      AVALC = coalesce(!!!select(., starts_with("AVAL") &
                                   (ends_with("C") |  ends_with("C.y"))))
    ) %>%
    select(-ends_with(".y"))

  print(add_data %>% select(USUBJID, PARAMCD, starts_with("AVA")))


  # add the dates for the derived parameters
  by_vars <- vars2chr(by_vars)
  dates <- subset_ds %>%
    group_by(!!!syms(by_vars)) %>%
    summarise(
      ASTDTM___ = min(ASTDTM, na.rm = TRUE),
      AENDTM___ = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
    )

  expo_data <- add_data %>%
    left_join(dates, by = by_vars) %>%
    mutate(
      ASTDTM = coalesce(as_iso_dttm(ASTDTM), as_iso_dttm(ASTDTM___)),
      AENDTM = coalesce(as_iso_dttm(AENDTM), as_iso_dttm(AENDTM___)),
      ASTDT = date(ASTDTM),
      AENDT = date(AENDTM),
      !!!set_values_to
    ) %>%
    select(-ends_with("___"))

  data <- bind_rows(dataset, expo_data)
}


