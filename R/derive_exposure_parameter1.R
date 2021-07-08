

derive_exposure_parameters <- function(dataset,
                                       by_vars,
                                       filter_rows = NULL,
                                       set_values_to = NULL,
                                       drop_values_from = NULL,
                                       ...) {
  assert_data_frame(dataset)
  assert_vars(drop_values_from, optional = TRUE)
  # need to add check that ASTDT(M)/AENDT(M) are present to derive the dates
  by_vars <- assert_vars(by_vars)
  filter_rows <- assert_filter_cond(enquo(filter_rows), optional = TRUE)
  if (!quo_is_null(filter_rows)) {
    subset_ds <- dataset %>%
      filter(!!filter_rows)
  } else {
    subset_ds <- dataset
  }

  if (!is.null(set_values_to)) {
    assert_that(is_quosures(set_values_to),
      msg = str_glue("`set_values_to` must be a `vars()` object, \\
                             not {friendly_type(typeof(set_values_to))}.")
    )
  }

  sources <- list(...)
  walk(sources, validate_exposure_var_spec)
  add_data <- vector("list", length(sources))

  for (i in seq_along(sources)) {

    input_param<-sources[[i]]$input_parameters
    new_param<-sources[[i]]$new_parameters

    print(input_param)
    print(new_param)

    add_data[[i]] <- subset_ds

    add_data_j <- vector("list", length(new_param))

    for (j in seq_along(new_param) ){
     add_data_j[[j]]<- add_data[[i]] %>%
          filter(PARAMCD == input_param[[j]])

      print(c(i,j))
      print(add_data_j[[j]] %>%
              select(USUBJID, PARAMCD))
      print(new_param[[j]])
      print(input_param[[j]])

      add_data_j[[j]] <- add_data_j[[j]] %>%
        derive_summary_records(
          by_vars = by_vars,
          filter_rows = PARAMCD %in% c(!!input_param[[j]]),
          fns = sources[[i]]$fns,
          set_values_to = vars(PARAMCD___ = !!new_param[[j]])  ,
          drop_values_from=drop_values_from
        ) %>%
        # keep only the new param
        filter(PARAMCD___ == !!new_param[[j]]) %>%
        mutate(PARAMCD=PARAMCD___)%>%
        select(-ends_with("___"))
        #select(-TEMP___)
      print(add_data_j[[j]] %>%
              select(USUBJID, PARAMCD))

  }

  add_data[[i]]<-bind_rows(add_data_j)

  print(head(add_data[[i]] ))

  }

  expo_data <- bind_rows(add_data)

  # add the dates for the derived parameters
  by_vars <- vars2chr(by_vars)
  dates <- subset_ds %>%
    group_by(!!!syms(by_vars)) %>%
    summarise(
      ASTDTM___ = min(ASTDTM, na.rm = TRUE),
      AENDTM___ = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE)
    )

  expo_data <- expo_data %>%
    left_join(dates, by = by_vars) %>%
    mutate(
      ASTDTM = coalesce(as_iso_dttm(ASTDTM), as_iso_dttm(ASTDTM___)),
      AENDTM = coalesce(as_iso_dttm(AENDTM), as_iso_dttm(AENDTM___)),
      ASTDT = date(ASTDTM),
      AENDT = date(AENDTM),
      #PARAMCD= PARAMCD___,
      !!!set_values_to) %>%
    select(-ends_with("___"))

  data <- bind_rows(dataset, expo_data)
}



exposure_var_spec <- function(new_parameters, input_parameters, fns) {
  out <- list(
    new_parameters = new_parameters,
    input_parameters = input_parameters,
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
  #assert_symbol(values$new_parameters)
  #assert_symbol(values$input_parameters)
  assert_character_vector(values$new_parameters)
  assert_character_vector(values$input_parameters)
  # check fns???
   obj
}




