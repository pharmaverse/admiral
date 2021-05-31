
#' Derive last dose date(-time)
#'
#' @param dataset Input AE dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' Defaults to NULL.
#' @param by_vars Variables to join by (type `vars`).
#' @param dose_start The dose start date variable.
#' @param dose_end The dose end date variable.
#' @param analysis_date The analysis date variable.
#' @param output_var The output variable.
#' @param output_datetime Logical. Should only date or date-time variable be returned?
#' Defaults to `TRUE` (i.e. date-time variable).
#' @param check_dates_only Logical.
#' An assumption that start and end dates of treatment match is checked.
#' By default (`FALSE`), the date as well as the time component is checked.
#' If set to `TRUE`, then only the date component of those variables is checked.
#'
#' @details All date (date-time) variables can be characters in standard ISO format or
#' of date / date-time class.
#' For ISO format, see [`impute_dtc`], parameter `dtc` for further details.
#'
#' The last dose date is derived as follows:
#' Firstly, the `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' Secondly, the datasets `dataset` and `dataset_ex` are joined using `by_vars`
#' and `AESEQ` variable.
#' Thirdly, the last dose date is derived:
#' the last dose date is the maximum date where `dose_end` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' Lastly, the last dose date is appended to the `dataset` and returned to the user.
#'
#' @return AE dataset with additional column `output_var`.
#'
#' @author Ondrej Slama
#'
#' @export
#'
#' @examples
#' data(ae); data(ex_single)
#' derive_last_dose(
#'   ae,
#'   ex_single,
#'   filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & stringr::str_detect(EXTRT, "PLACEBO"))) &
#'     nchar(as.character(EXENDTC)) >= 10,
#'   dose_start = EXSTDTC,
#'   dose_end = EXENDTC,
#'   analysis_date = AESTDTC,
#'   output_var = LDOSEDTM,
#'   output_datetime = TRUE,
#'   check_dates_only = FALSE
#' )
#'
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = NULL,
                             by_vars = vars(STUDYID, USUBJID),
                             dose_start,
                             dose_end,
                             analysis_date,
                             output_var,
                             output_datetime = TRUE,
                             check_dates_only = FALSE) {

  assert_that(
    !dplyr::is_grouped_df(dataset),
    !dplyr::is_grouped_df(dataset_ex),
    is.list(by_vars),
    rlang::is_scalar_logical(output_datetime),
    rlang::is_scalar_logical(check_dates_only)
  )

  # apply filtering condition
  filter_ex <- enquo(filter_ex)
  if (!is.null(rlang::quo_get_expr(filter_ex))) {
    dataset_ex <- filter(dataset_ex, !!filter_ex)
  }

  dose_start <- enquo(dose_start)
  dose_end <- enquo(dose_end)
  analysis_date <- enquo(analysis_date)
  output_var <- enquo(output_var)
  by_vars_str <- vapply(by_vars,
                        function(x) as_string(rlang::quo_get_expr(x)),
                        character(1),
                        USE.NAMES = F)

  # check variables existence
  assert_has_variables(
    dataset,
    c(by_vars_str,
      "AESEQ",
      as_string(rlang::quo_get_expr(analysis_date))
    )
  )

  assert_has_variables(
    dataset_ex,
    c(by_vars_str,
      as_string(rlang::quo_get_expr(dose_start)),
      as_string(rlang::quo_get_expr(dose_end)))
  )

  # assumption for last dose derivation: start and end dates (datetimes) need to match
  if (check_dates_only) {
    check_cond <- dplyr::summarise(dataset_ex,
                                   all_equal = all(as.Date(!!dose_start) == as.Date(!!dose_end)))
  } else {
    check_cond <- dplyr::summarise(dataset_ex,
                                   all_equal = all(!!dose_start == !!dose_end))
  }
  if (!check_cond$all_equal) {
    stop(paste(
      "Not all values of", as_string(rlang::quo_get_expr(dose_start)),
      "are equal to", as_string(rlang::quo_get_expr(dose_end))
    ))
  }

  # select only a subset of columns
  dataset_ex <- select(dataset_ex, !!!by_vars, !!dose_end)

  # calculate last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    dplyr::mutate_at(dplyr::vars(!!dose_end, !!analysis_date),
                     ~ `if`(is_date(.), convert_dtm_to_dtc(.), .)) %>%
    group_by(!!!by_vars, .data$AESEQ) %>%
    mutate(
      tmp_exendtc = impute_dtc(dtc = !!dose_end,
                               date_imputation = NULL,
                               time_imputation = "00:00:00") %>%
        convert_dtc_to_dtm(),
      tmp_aestdtc = impute_dtc(dtc = !!analysis_date,
                               date_imputation = NULL,
                               time_imputation = "23:59:59") %>%
        convert_dtc_to_dtm()) %>%
    dplyr::summarise(!!output_var := compute_ldosedtm(exendtc = .data$tmp_exendtc,
                                                      aestdtc = .data$tmp_aestdtc)) %>%
    ungroup()

  # return either date or date-time variable
  if (!output_datetime) {
    res <- mutate(res, !!output_var := as.Date(!!output_var))
  }

  # return dataset with additional column
  left_join(dataset,
            dplyr::distinct(res, !!!by_vars, .data$AESEQ, !!output_var),
            by = c(by_vars_str, "AESEQ"))
}

#' Helper function to calculate last dose
#'
#' @param exendtc dose end date
#' @param aestdtc analysis date
#'
#' @return date-time vector
compute_ldosedtm <- function(exendtc, aestdtc) {
  if (any(!is.na(exendtc) & !is.na(aestdtc)) && any(exendtc <= aestdtc)) {
    max(exendtc[exendtc <= aestdtc])
  } else {
    as.POSIXct(NA)
  }
}

#' Helper function to convert date (or date-time) objects to characters of dtc format
#' (-DTC type of variable)
#'
#' @param dtm date or date-time
#'
#' @return character
#'
#' @examples
#' admiral:::convert_dtm_to_dtc(as.POSIXct(Sys.time()))
#' admiral:::convert_dtm_to_dtc(as.Date(Sys.time()))
convert_dtm_to_dtc <- function(dtm) {
  stopifnot(lubridate::is.instant(dtm))
  format(dtm, "%Y-%m-%dT%H:%M:%S")
}
