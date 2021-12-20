#' Derive Single Dose
#'
#' Derives Single Dose from aggregate dose
#'
#' @param dataset Input dataset
#'
#'   The columns specified by `by_vars`, `dose_freq`, `start_date` and the `end_date`
#'   parameters are expected.
#'
#' @param by_vars Variables to group by (created by `dplyr::vars`).
#'
#'   These variables are used to identify unique rows for each aggregate dose.
#'
#'   Example: vars(STUDYID, USUBJID, EXTRT)
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `ASTDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `AENDT`
#'
#' @param dose_freq The dose frequency
#'
#'   The aggregate dosing frequency used for multiple doses in a row.
#'
#'   Default: EXDOSFRQ
#'
#'   Permitted Values: 'QxD', 'QxW', 'EVERY x DAYS', 'EVERY x WEEKS'
#'   where x is a positive integer
#'
#' @details Each aggregate dose row is split into multiple rows which each represent
#' a single dose.
#'
#' @author Michael Thorpe
#'
#' @return The input dataset with a single dose per row.
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
#'   "P01", "Q2D", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-07"),
#'   "P01", "Q3D", lubridate::ymd("2021-01-08"), lubridate::ymd("2021-01-15"),
#'   "P01", "EVERY 2 WEEKS", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-29")
#' )
#'
#' derive_vars_single_dose(data, by_vars = vars(USUBJID))
derive_vars_single_dose <- function(dataset,
                                    by_vars,
                                    dose_freq = EXDOSFRQ,
                                    start_date = ASTDT,
                                    end_date = AENDT) {

  col_names <- colnames(dataset)
  by_vars <- assert_vars(by_vars)
  dose_freq <- assert_symbol(enquo(dose_freq))
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, required_vars = quo_c(by_vars, dose_freq, start_date, end_date))

  pattern = c("D", "W")

  dataset <- dataset %>%
    mutate(temp_freq = as.numeric(str_extract(!!dose_freq, "[:digit:]")),
           temp_unit = str_extract(!!dose_freq, paste(pattern, collapse = "|")),
           time_gap = as.numeric(!!end_date - !!start_date),
           new_dose_no = case_when(temp_unit == "D" ~ 1,
                                   temp_unit == "W" ~ 7),
           dose_freq = temp_freq*new_dose_no,
           num_of_doses = floor(time_gap/(dose_freq)) + 1
    )


  dataset <- dataset[rep(row.names(dataset), dataset$num_of_doses),]

  dataset <- dataset %>%
    group_by(!!!by_vars, !!dose_freq, !!start_date, !!end_date) %>%
    mutate(temp_multiplier = (row_number() - 1)) %>%
    ungroup() %>%
    mutate(day_difference = days(temp_multiplier*dose_freq))

  dataset <- dataset %>%
    mutate(!!dose_freq := "ONCE",
           !!end_date := !!start_date + day_difference,
           !!start_date := !!start_date + day_difference
    ) %>%
    select(all_of(col_names))

  return(dataset)
}

