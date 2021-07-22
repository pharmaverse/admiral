

derive_exposure_params <- function(dataset,
                                   by_vars,
                                   filter_rows = NULL,
                                   set_values_to = NULL,
                                   drop_values_from = NULL,
                                   ...) {
  assert_data_frame(dataset, required_vars = vars(ASTDTM, AENDTM))
  assert_vars(drop_values_from, optional = TRUE)
  by_vars <- assert_vars(by_vars)
  filter_rows <- assert_filter_cond(enquo(filter_rows), optional = TRUE)
  if (!is.null(set_values_to)) {
    assert_that(is_quosures(set_values_to),
      msg = str_glue("`set_values_to` must be a `vars()` object, \\
                             not {friendly_type(typeof(set_values_to))}.")
    )
  }

  subset_ds <- dataset %>%
    filter_if(filter_rows)

  sources <- list(...)
  s1 <- sources[sapply(sources, function(x) class(x)[1] == "exposure_var_spec")]

  walk(s1, validate_exposure_var_spec)
  add_data <- vector("list", length(s1))
  for (i in seq_along(s1)) {
    add_data[[i]] <- subset_ds %>%
      filter(PARAMCD == s1[[i]]$input_parameter) %>%
      derive_summary_records(
        by_vars = by_vars,
        fns = s1[[i]]$fns,
        set_values_to = vars(PARAMCD___ = !!s1[[i]]$new_parameter),
        drop_values_from = drop_values_from
      ) %>%
      filter(PARAMCD___ == !!s1[[i]]$new_parameter) %>%
      mutate(PARAMCD = PARAMCD___) %>%
      select(-ends_with("___"))
  }
  expo_data1 <- bind_rows(add_data) %>%
    # not sure why the derive_summary_records() fns render a AVAL.x and AVAL.y...
    # AVAL.x retain the value of the input param, AVAL.y has the correctresult
    select(-ends_with(".x")) %>%
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


  # add the dates for the derived parameters
  by_vars <- vars2chr(by_vars)
  dates <- subset_ds %>%
    group_by(!!!syms(by_vars)) %>%
    summarise(
      ASTDTM___ = min(ASTDTM, na.rm = TRUE),
      AENDTM___ = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
    )

  expo_data <- expo_data1 %>%
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



exposure_var_spec <- function(new_parameter, input_parameter, fns) {
  out <- list(
    new_parameter = new_parameter,
    input_parameter = input_parameter,
    fns = fns
  )
  class(out) <- c("exposure_var_spec", "list")
  validate_exposure_var_spec(out)
}

#' Validate an object is indeed a `_source` object
#'
#' @param obj An object to be validated.
#'
#' @author Samia Kabi
#'
#' @noRd
#'
#' @return The original object.
validate_exposure_var_spec <- function(obj) {
  assert_that(inherits(obj, "exposure_var_spec"))
  values <- unclass(obj)
  assert_character_scalar(values$new_parameter)
  assert_character_scalar(values$input_parameter)
  # check fns???
  obj
}

