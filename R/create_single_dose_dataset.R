#' Create dataset of single doses
#'
#' Derives dataset of single dose from aggregate dose information. This may be
#' necessary when e.g. calculating last dose before an adverse event in `ADAE`
#' or deriving a total dose parameter in `ADEX` when `EXDOSFRQ != ONCE`.
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("dose_freq", "start_date", "end_date"))`
#'
#' @param dose_freq The dose frequency
#'
#'   The aggregate dosing frequency used for multiple doses in a row.
#'
#'   Permitted Values: defined by lookup table.
#'
#' @param start_date The start date
#'
#'   A date object is expected. This object cannot contain `NA` values.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param start_datetime The start date-time
#'
#'   A date-time object is expected. This object cannot contain `NA` values.
#'
#'   Refer to `derive_vars_dtm()` to impute and derive a date-time from a date
#'   character vector to a date object.
#'
#'   If the input dataset contains frequencies which refer to `DOSE_WINDOW`
#'   equals `"HOUR"` or `"MINUTE"`, the parameter must be specified.
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected. This object cannot contain `NA` values.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param end_datetime The end date-time
#'
#'   A date-time object is expected. This object cannot contain `NA` values.
#'
#'   Refer to `derive_vars_dtm()` to impute and derive a date-time from a date
#'   character vector to a date object.
#'
#'   If the input dataset contains frequencies which refer to `DOSE_WINDOW`
#'   equals `"HOUR"` or `"MINUTE"`, the parameter must be specified.
#'
#' @param lookup_table The dose frequency value lookup table
#'
#'   The table used to look up `dose_freq` values and determine the appropriate
#'   multiplier to be used for row generation. If a lookup table other than the
#'   default is used, it must have columns `DOSE_WINDOW`, `DOSE_COUNT`, and
#'   `CONVERSION_FACTOR`. The default table `dose_freq_lookup` is described in
#'   detail [here][dose_freq_lookup].
#'
#'   Permitted Values for `DOSE_WINDOW`: `"MINUTE"`, `"HOUR"`, `"DAY"`,
#'   `"WEEK"`, `"MONTH"`, `"YEAR"`
#'
#' @param lookup_column The dose frequency value column in the lookup table
#'
#'   The column of `lookup_table`.
#'
#' @param nominal_time The nominal relative time from first dose (`NFRLT`)
#'
#'   Used for PK analysis, this will be in hours and should be 0 for
#'   the first dose.  It can be derived as `(VISITDY - 1) * 24` for example.
#'   This will be expanded as the single dose dataset is created.  For example
#'   an `EXDOFRQ` of `"QD"` will result in the nominal_time being incremented by
#'   24 hours for each expanded record.
#'
#'   The value can be NULL if not needed.
#'
#' @param keep_source_vars List of variables to be retained from source dataset
#'
#'   This parameter can be specified if additional information is required in
#'   the output dataset. For example `EXTRT` for studies with more than one
#'   drug.
#'
#' @details Each aggregate dose row is split into multiple rows which each
#'   represent a single dose.The number of completed dose periods between
#'   `start_date` or `start_datetime` and `end_date` or `end_datetime` is
#'   calculated with `compute_duration` and multiplied by `DOSE_COUNT`.
#'   For `DOSE_WINDOW` values of `"WEEK"`, `"MONTH"`, and `"YEAR"`,
#'   `CONVERSION_FACTOR` is used to convert into days the time object
#'   to be added to `start_date`.
#'
#'   Observations with dose frequency `"ONCE"` are copied to the output dataset
#'   unchanged.
#'
#'
#' @family create_aux
#' @keywords create_aux
#'
#' @return The input dataset with a single dose per row.
#'
#' @export
#'
#' @examples
#' # Example with default lookup
#'
#' library(lubridate)
#' library(stringr)
#' library(tibble)
#' library(dplyr)
#'
#' data <- tribble(
#'   ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
#'   "STUDY01", "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
#'   ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
#'   "STUDY01", "P01", "Q3D", ymd("2021-01-08"), ymd_hms("2021-01-08 12:00:00"),
#'   ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
#'   "STUDY01", "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
#'   ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
#' )
#'
#' create_single_dose_dataset(data)
#'
#' # Example with custom lookup
#'
#' custom_lookup <- tribble(
#'   ~Value,   ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
#'   "Q30MIN", (1 / 30),    "MINUTE",                      1,
#'   "Q90MIN", (1 / 90),    "MINUTE",                      1
#' )
#'
#' data <- tribble(
#'   ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
#'   "STUDY01", "P01", "Q30MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
#'   ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
#'   "STUDY02", "P02", "Q90MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
#'   ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00")
#' )
#'
#' create_single_dose_dataset(data,
#'   lookup_table = custom_lookup,
#'   lookup_column = Value,
#'   start_datetime = ASTDTM,
#'   end_datetime = AENDTM
#' )
#' # Example with nominal time
#'
#' data <- tribble(
#'   ~STUDYID, ~USUBJID, ~EXDOSFRQ, ~NFRLT, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
#'   "STUDY01", "P01", "BID", 0, ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
#'   ymd("2021-01-07"), ymd_hms("2021-01-07 20:00:00"),
#'   "STUDY01", "P01", "BID", 168, ymd("2021-01-08"), ymd_hms("2021-01-08 08:00:00"),
#'   ymd("2021-01-14"), ymd_hms("2021-01-14 20:00:00"),
#'   "STUDY01", "P01", "BID", 336, ymd("2021-01-15"), ymd_hms("2021-01-15 08:00:00"),
#'   ymd("2021-01-29"), ymd_hms("2021-01-29 20:00:00")
#' )
#'
#' create_single_dose_dataset(data,
#'   dose_freq = EXDOSFRQ,
#'   start_date = ASTDT,
#'   start_datetime = ASTDTM,
#'   end_date = AENDT,
#'   end_datetime = AENDTM,
#'   lookup_table = dose_freq_lookup,
#'   lookup_column = CDISC_VALUE,
#'   nominal_time = NFRLT,
#'   keep_source_vars = exprs(
#'     USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, NFRLT
#'   )
#' )
#'
#' # Example - derive a single dose dataset with imputations
#'
#' # For either single drug administration records, or multiple drug administration
#' # records covering a range of dates, fill-in of missing treatment end datetime
#' # `EXENDTC` by substitution with an acceptable alternate, for example date of
#' # death, date of datacut may be required. This example shows the
#' # maximum possible number of single dose records to be derived. The example
#' # requires the date of datacut `DCUTDT` to be specified correctly, or
#' # if not appropriate to use `DCUTDT` as missing treatment end data and missing
#' # treatment end datetime could set equal to treatment start date and treatment
#' # start datetime. ADSL variables `DTHDT` and `DCUTDT` are preferred for
#' # imputation use.
#' #
#' # All available trial treatments are included, allowing multiple different
#' # last dose variables to be created in for example `use_ad_template("ADAE")`
#' # if required.
#'
#' adsl <- tribble(
#'   ~STUDYID, ~USUBJID, ~DTHDT,
#'   "01", "1211", ymd("2013-01-14"),
#'   "01", "1083", ymd("2013-08-02"),
#'   "01", "1445", ymd("2014-11-01"),
#'   "01", "1015", NA,
#'   "01", "1023", NA
#' )
#'
#' ex <- tribble(
#'   ~STUDYID, ~USUBJID, ~EXSEQ, ~EXTRT, ~EXDOSE, ~EXDOSU, ~EXDOSFRQ, ~EXSTDTC, ~EXENDTC,
#'   "01", "1015", 1, "PLAC", 0, "mg", "QD", "2014-01-02", "2014-01-16",
#'   "01", "1015", 2, "PLAC", 0, "mg", "QD", "2014-06-17", "2014-06-18",
#'   "01", "1015", 3, "PLAC", 0, "mg", "QD", "2014-06-19", NA_character_,
#'   "01", "1023", 1, "PLAC", 0, "mg", "QD", "2012-08-05", "2012-08-27",
#'   "01", "1023", 2, "PLAC", 0, "mg", "QD", "2012-08-28", "2012-09-01",
#'   "01", "1211", 1, "XANO", 54, "mg", "QD", "2012-11-15", "2012-11-28",
#'   "01", "1211", 2, "XANO", 54, "mg", "QD", "2012-11-29", NA_character_,
#'   "01", "1445", 1, "PLAC", 0, "mg", "QD", "2014-05-11", "2014-05-25",
#'   "01", "1445", 2, "PLAC", 0, "mg", "QD", "2014-05-26", "2014-11-01",
#'   "01", "1083", 1, "PLAC", 0, "mg", "QD", "2013-07-22", "2013-08-01"
#' )
#'
#' adsl_death <- adsl %>%
#'   mutate(
#'     DTHDTM = convert_date_to_dtm(DTHDT),
#'     # Remove `DCUT` setup line below if ADSL `DCUTDT` is populated.
#'     DCUTDT = convert_dtc_to_dt("2015-03-06"), # Example only, enter date.
#'     DCUTDTM = convert_date_to_dtm(DCUTDT)
#'   )
#'
#' # Select valid dose records, non-missing `EXSTDTC` and `EXDOSE`.
#' ex_mod <- ex %>%
#'   filter(!is.na(EXSTDTC) & !is.na(EXDOSE)) %>%
#'   derive_vars_merged(adsl_death, by_vars = get_admiral_option("subject_keys")) %>%
#'   # Example, set up missing `EXDOSFRQ` as QD daily dosing regime.
#'   # Replace with study dosing regime per trial treatment.
#'   mutate(EXDOSFRQ = if_else(is.na(EXDOSFRQ), "QD", EXDOSFRQ)) %>%
#'   # Create EXxxDTM variables and replace missing `EXENDTM`.
#'   derive_vars_dtm(
#'     dtc = EXSTDTC,
#'     new_vars_prefix = "EXST",
#'     date_imputation = "first",
#'     time_imputation = "first",
#'     flag_imputation = "none",
#'   ) %>%
#'   derive_vars_dtm_to_dt(exprs(EXSTDTM)) %>%
#'   derive_vars_dtm(
#'     dtc = EXENDTC,
#'     new_vars_prefix = "EXEN",
#'     # Maximum imputed treatment end date must not be not greater than
#'     # date of death or after the datacut date.
#'     max_dates = exprs(DTHDTM, DCUTDTM),
#'     date_imputation = "last",
#'     time_imputation = "last",
#'     flag_imputation = "none",
#'     highest_imputation = "Y",
#'   ) %>%
#'   derive_vars_dtm_to_dt(exprs(EXENDTM)) %>%
#'   # Select only unique values.
#'   # Removes duplicated records before final step.
#'   distinct(
#'     STUDYID, USUBJID, EXTRT, EXDOSE, EXDOSFRQ, DCUTDT, DTHDT, EXSTDT,
#'     EXSTDTM, EXENDT, EXENDTM, EXSTDTC, EXENDTC
#'   )
#'
#' create_single_dose_dataset(
#'   ex_mod,
#'   start_date = EXSTDT,
#'   start_datetime = EXSTDTM,
#'   end_date = EXENDT,
#'   end_datetime = EXENDTM,
#'   keep_source_vars = exprs(
#'     STUDYID, USUBJID, EXTRT, EXDOSE, EXDOSFRQ,
#'     DCUTDT, EXSTDT, EXSTDTM, EXENDT, EXENDTM, EXSTDTC, EXENDTC
#'   )
#' )
create_single_dose_dataset <- function(dataset,
                                       dose_freq = EXDOSFRQ,
                                       start_date = ASTDT,
                                       start_datetime = NULL,
                                       end_date = AENDT,
                                       end_datetime = NULL,
                                       lookup_table = dose_freq_lookup,
                                       lookup_column = CDISC_VALUE,
                                       nominal_time = NULL,
                                       keep_source_vars = expr_c(
                                         get_admiral_option("subject_keys"), dose_freq,
                                         start_date, start_datetime,
                                         end_date, end_datetime
                                       )) {
  dose_freq <- assert_symbol(enexpr(dose_freq))
  lookup_column <- assert_symbol(enexpr(lookup_column))
  start_date <- assert_symbol(enexpr(start_date))
  start_datetime <- assert_symbol(enexpr(start_datetime), optional = TRUE)
  end_date <- assert_symbol(enexpr(end_date))
  end_datetime <- assert_symbol(enexpr(end_datetime), optional = TRUE)
  nominal_time <- assert_symbol(enexpr(nominal_time), optional = TRUE)
  assert_data_frame(dataset, required_vars = expr_c(dose_freq, start_date, end_date))
  assert_data_frame(
    lookup_table,
    required_vars = exprs(!!lookup_column, DOSE_WINDOW, DOSE_COUNT, CONVERSION_FACTOR)
  )
  assert_data_frame(dataset, required_vars = keep_source_vars)
  col_names <- colnames(dataset)

  # Checking that the dates specified follow the ADaM naming convention of ending in DT
  start_datec <- as_string(as_name(start_date))
  start_date_chk <- str_locate_all(start_datec, "DT")
  start_date_chk_pos <- as.vector(start_date_chk[[1]])

  if (str_length(start_datec) != start_date_chk_pos[-1]) {
    err_msg <- c(
      "The argument {.arg start_date} is expected to have a name  ending with {.val ---DT}.",
      "Please update as it does not follow the expected naming convention."
    )
    cli_abort(err_msg)
  }

  end_datec <- as_string(as_name(end_date))
  end_date_chk <- str_locate_all(end_datec, "DT")
  end_date_chk_pos <- as.vector(end_date_chk[[1]])

  if (str_length(end_datec) != end_date_chk_pos[-1]) {
    err_msg <- c(
      "The argument {.arg end_date} is expected to have a name like {.val xxxDT}.",
      "Please check as it does not follow the expected naming convention."
    )
    cli_abort(err_msg)
  }

  # Set up lookup table to be joined to dataset

  lookup <- lookup_table %>%
    rename(!!dose_freq := !!lookup_column)

  # Observations with frequency ONCE are copied unchanged to the output dataset
  data_once <- filter(dataset, !!dose_freq == "ONCE")

  data_not_once <- filter(dataset, !!dose_freq != "ONCE")

  # Check that NAs do not appear in start_date or start_datetime or end_date or end_datetime columns
  condition <- paste0("is.na(", as_name(start_date), ") | is.na(", as_name(end_date), ")")
  if (!is.null(start_datetime)) {
    condition <- paste(condition, "| is.na(", as_name(start_datetime), ")")
  }
  if (!is.null(end_datetime)) {
    condition <- paste(condition, "| is.na(", as_name(end_datetime), ")")
  }
  na_check <- data_not_once %>%
    filter(!!parse_expr(condition)) %>%
    select(!!start_date, !!end_date, !!start_datetime, !!end_datetime)

  if (nrow(na_check) > 0) {
    na_columns <- paste0(colnames(na_check)[colSums(is.na(na_check)) > 0], collapse = ", ")
    err_msg <- c(
      "The arguments {.arg start_date} or {.arg start_datetime} and
       {.arg end_date} or {.arg end_datetime} cannot contain {.val {NA}} values.",
      "i" = "Please check {.val {na_columns}} for {.val {NA}} values."
    )

    cli_abort(err_msg)
  }

  # Check nominal time is NULL or numeric

  if (!is.null(nominal_time)) {
    check_nominal_time <-
      mutate(
        data_not_once,
        check_nom := !!nominal_time
      )

    assert_numeric_vector(check_nominal_time$check_nom)
  }

  # Check values of lookup vs. data and return error if values are not covered

  value_check <- data_not_once %>%
    select(!!dose_freq) %>%
    anti_join(lookup, by = as.character(dose_freq)) %>%
    unique()

  if (nrow(value_check) > 0) {
    values_not_found <- value_check %>% select(!!dose_freq)

    err_msg <-
      c("The following values of {.val {as.character(dose_freq)}} in
         {.arg dataset} do not appear in
         {.arg lookup_table}: {.val {values_not_found}}")
    cli_abort(err_msg)
  }

  # Check lookup_table does not contain duplicates
  signal_duplicate_records(
    dataset = lookup_table,
    by_vars = exprs(!!lookup_column)
  )

  # Use compute_duration to determine the number of completed dose periods

  if (is.null(start_datetime)) {
    min_hour_cases <- exprs(FALSE ~ 0)
  } else {
    min_hour_cases <- exprs(
      DOSE_WINDOW == "MINUTE" ~ compute_duration(!!start_datetime, !!end_datetime,
        in_unit = "minutes", out_unit = "minutes"
      ),
      DOSE_WINDOW == "HOUR" ~ compute_duration(!!start_datetime, !!end_datetime,
        in_unit = "hours", out_unit = "hours"
      )
    )
  }
  data_not_once <- left_join(
    data_not_once,
    lookup,
    by = as.character(dose_freq)
  )

  if (any(data_not_once$DOSE_WINDOW %in% c("MINUTE", "HOUR")) &&
    (is.null(start_datetime) || is.null(end_datetime))) {
    cli_abort(
      "There are dose frequencies more frequent than once a day, thus
       arguments {.arg start_datetime} and {.arg end_datetime} must be specified."
    )
  }

  data_not_once <- data_not_once %>%
    mutate(dose_periods = case_when(
      !!!min_hour_cases,
      DOSE_WINDOW == "DAY" ~ compute_duration(!!start_date, !!end_date, out_unit = "days"),
      DOSE_WINDOW == "WEEK" ~ compute_duration(!!start_date, !!end_date, out_unit = "weeks"),
      DOSE_WINDOW == "MONTH" ~ compute_duration(!!start_date, !!end_date, out_unit = "months"),
      DOSE_WINDOW == "YEAR" ~ compute_duration(!!start_date, !!end_date, out_unit = "years")
    )) %>%
    mutate(dose_count = ceiling(dose_periods * DOSE_COUNT)) %>%
    derive_var_obs_number(new_var = grpseq)

  # Generate a row for each completed dose

  data_not_once <- data_not_once[rep(row.names(data_not_once), data_not_once$dose_count), ]

  # Determine amount of days to adjust start_date or start_datetime and end_date or end_datetime

  data_not_once <- data_not_once %>%
    group_by(grpseq, !!dose_freq, !!start_date, !!end_date) %>%
    mutate(time_increment = (row_number() - 1) / (DOSE_COUNT)) %>%
    ungroup() %>%
    mutate(
      time_differential = case_when(
        DOSE_WINDOW == "MINUTE" ~ minutes(floor(time_increment)),
        DOSE_WINDOW == "HOUR" ~ hours(floor(time_increment)),
        DOSE_WINDOW %in% c("DAY", "WEEK", "MONTH", "YEAR") ~
          days(floor(time_increment / CONVERSION_FACTOR))
      ),
      time_differential_dt = case_when(
        DOSE_WINDOW == "MINUTE" ~ days(floor(time_increment / 1440)),
        DOSE_WINDOW == "HOUR" ~ days(floor(time_increment / 24)),
        DOSE_WINDOW %in% c("DAY", "WEEK", "MONTH", "YEAR") ~
          days(floor(time_increment / CONVERSION_FACTOR))
      )
    )

  # Adjust start_date and end_date, drop calculation columns, make sure nothing
  # later than end_date shows up in output

  data_not_once <- data_not_once %>%
    mutate(
      !!dose_freq := "ONCE",
      !!start_date := !!start_date + time_differential_dt
    )
  if (!is.null(start_datetime)) {
    data_not_once <-
      mutate(
        data_not_once,
        !!start_datetime := !!start_datetime + time_differential,
      )
  }
  if (!is.null(nominal_time)) {
    data_not_once <-
      mutate(
        data_not_once,
        !!nominal_time := !!nominal_time + as.numeric(time_differential) / 3600
      )
  }

  data_not_once <- data_not_once %>%
    filter(!(!!start_date > !!end_date)) %>%
    mutate(
      !!end_date := !!start_date
    )
  if (!is.null(end_datetime)) {
    data_not_once <-
      mutate(
        data_not_once,
        !!end_datetime := case_when(
          DOSE_WINDOW %in% c("MINUTE", "HOUR") ~ !!start_datetime,
          DOSE_WINDOW %in% c("DAY", "WEEK", "MONTH", "YEAR") ~
            ymd_hms(paste0(!!start_date, " ", format(!!end_datetime, format = "%H:%M:%S")))
        )
      )
  }
  data_not_once <- select(data_not_once, all_of(col_names))

  # Stitch back together

  bind_rows(data_once, data_not_once) %>%
    select(!!!keep_source_vars)
}
