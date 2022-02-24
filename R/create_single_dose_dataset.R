#' Pre-Defined Dose Frequencies
#'
#' These pre-defined dose frequencies are sourced from
#' [CDISC](https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.pdf). The
#' number of rows to generate using `create_single_dose_dataset()` arguments
#' `start_date` and `end_date` is derived from `DFFactor` and `DFUnit` with
#' appropriate functions from `lubridate`.
#'
#' @details
#' To see the entire table in the console, run `print(dose_freq_lookup)`.
#'
#' @seealso [create_single_dose_dataset()]
#'
#' @export
#'
#' @keywords dose_freq
#'
#' @rdname dose_freq_lookup
#'
dose_freq_lookup <- tibble::tribble(
~NCICode, ~CDISCSubmissionValue
"C64526","1 TIME PER WEEK",
"C139179","10 DAYS PER MONTH",
"C176288","2 TIMES PER CYCLE",
"C64497","2 TIMES PER WEEK",
"C98861","2 TIMES PER YEAR",
"C176289","3 TIMES PER CYCLE",
"C98859","3 TIMES PER MONTH",
"C64528","3 TIMES PER WEEK",
"C98860","3 TIMES PER YEAR",
"C98852","4 TIMES PER MONTH",
"C64531","4 TIMES PER WEEK",
"C98853","4 TIMES PER YEAR",
"C98849","5 TIMES PER DAY",
"C98850","5 TIMES PER MONTH",
"C85552","5 TIMES PER WEEK",
"C98851","5 TIMES PER YEAR",
"C98855","6 TIMES PER DAY",
"C98856","6 TIMES PER MONTH",
"C98857","6 TIMES PER WEEK",
"C98858","6 TIMES PER YEAR",
"C139180","7 TIMES PER DAY",
"C98854","7 TIMES PER WEEK",
"C139181","8 TIMES PER DAY",
"C139182","9 TIMES PER DAY",
"C64496","BID",
"C71129","BIM",
"C161332","EVERY 12 WEEKS",
"C161336","EVERY 16 WEEKS",
"C71127","EVERY 2 WEEKS",
"C64535","EVERY 3 WEEKS",
"C161333","EVERY 3 YEARS",
"C64529","EVERY 4 WEEKS",
"C103390","EVERY 5 WEEKS",
"C161334","EVERY 5 YEARS",
"C89788","EVERY 6 WEEKS"
) %>%
mutate(DoseCount = case_when(
  str_detect(CDISCSubmissionValue, "PER") ~ as.numeric(str_remove_all(CDISCSubmissionValue,"[\\D]")),
  str_detect(CDISCSubmissionValue, "EVERY") ~ 1/as.numeric(str_remove_all(CDISCSubmissionValue,"[\\D]")),
  str_detect(CDISCSubmissionValue, "BID|BIM") ~ 2),
  DoseWindow = case_when(
    str_detect(CDISCSubmissionValue, "EVERY|PER") ~ str_remove_all(sub(".* (\\w+)$", "\\1", x),"S"),
    CDISCSubmissionValue=="BID" ~ str_remove_all(sub(".* (\\w+)$", "\\1", x),"S"),

  )
       )

#' Create dataset of single doses
#'
#' Derives dataset of single dose from aggregate dose information
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
#' create_single_dose_dataset(data, by_vars = vars(USUBJID))
create_single_dose_dataset <- function(dataset,
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
  lookup <- dose_freq_lookup %>%
    rename(!!dose_freq := CDISCSubmissionValue)

  dataset <- dataset %>%
    left_join(lookup, by = as.character(dose_freq)[2]) %>%
    mutate(dose_windows = case_when(
      DFUnit == "Day" ~ compute_duration(!!start_date, !!end_date, out_unit = "days"),
      DFUnit == "Week" ~ compute_duration(!!start_date, !!end_date, out_unit = "weeks"),
      DFUnit == "Month" ~ compute_duration(!!start_date, !!end_date, out_unit = "months"),
      DFUnit == "Year" ~ compute_duration(!!start_date, !!end_date, out_unit = "years")
    )) %>%
    mutate(dose_count = floor(dose_windows*DFFactor) + (floor(dose_windows*DFFactor)==0))


  dataset <- dataset[rep(row.names(dataset), dataset$dose_count), ]

  dataset <- dataset %>%
    group_by(!!!by_vars, !!dose_freq, !!start_date, !!end_date) %>%
    mutate(temp_dose_multiplier = (row_number() - 1)/DFFactor) %>%
    ungroup() %>%
    mutate(time_differential = case_when(
      DFUnit == "Day" ~ days(temp_dose_multiplier),
      DFUnit == "Week" ~ weeks(temp_dose_multiplier),
      DFUnit == "Month" ~ months(temp_dose_multiplier),
      DFUnit == "Year" ~ years(temp_dose_multiplier)
      ))

  dataset <- dataset %>%
    mutate(!!dose_freq := "ONCE",
           !!end_date := !!start_date + time_differential,
           !!start_date := !!start_date + time_differential
    ) %>%
    select(all_of(!!!vars(col_names)))
  return(dataset)
}
