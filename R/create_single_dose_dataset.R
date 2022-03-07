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
  ~NCICode, ~CDISCSubmissionValue,
  "C64526", "1 TIME PER WEEK",
  "C139179", "10 DAYS PER MONTH",
#  "C176288", "2 TIMES PER CYCLE",
  "C64497", "2 TIMES PER WEEK",
  "C98861", "2 TIMES PER YEAR",
#  "C176289", "3 TIMES PER CYCLE",
  "C98859", "3 TIMES PER MONTH",
  "C64528", "3 TIMES PER WEEK",
  "C98860", "3 TIMES PER YEAR",
  "C98852", "4 TIMES PER MONTH",
  "C64531", "4 TIMES PER WEEK",
  "C98853", "4 TIMES PER YEAR",
  "C98849", "5 TIMES PER DAY",
  "C98850", "5 TIMES PER MONTH",
  "C85552", "5 TIMES PER WEEK",
  "C98851", "5 TIMES PER YEAR",
  "C98855", "6 TIMES PER DAY",
  "C98856", "6 TIMES PER MONTH",
  "C98857", "6 TIMES PER WEEK",
  "C98858", "6 TIMES PER YEAR",
  "C139180", "7 TIMES PER DAY",
  "C98854", "7 TIMES PER WEEK",
  "C139181", "8 TIMES PER DAY",
  "C139182", "9 TIMES PER DAY",
  "C64496", "BID",
  "C71129", "BIM",
  "C161332", "EVERY 12 WEEKS",
  "C161336", "EVERY 16 WEEKS",
  "C71127", "EVERY 2 WEEKS",
  "C64535", "EVERY 3 WEEKS",
  "C161333", "EVERY 3 YEARS",
  "C64529", "EVERY 4 WEEKS",
  "C103390", "EVERY 5 WEEKS",
  "C161334", "EVERY 5 YEARS",
  "C89788", "EVERY 6 WEEKS",
  "C116149",  "EVERY 7 WEEKS",
  "C103389",  "EVERY 8 WEEKS",
  "C154484",  "EVERY AFTERNOON",
  "C160957",  "EVERY EVENING",
  "C67069",  "EVERY WEEK",
  "C74924",  "PA",
  #"C64499",  "PRN",
  # "C64500",  "Q10H",
  # "C64501",  "Q11H",
  # "C64502",  "Q12H",
  # "C64503",  "Q13H",
  # "C64504",  "Q14H",
  # "C64505",  "Q15H",
  # "C64506",  "Q16H",
  # "C64507",  "Q17H",
  # "C64508",  "Q18H",
  # "C64509",  "Q19H",
  # "C64511",  "Q20H",
  # "C64512",  "Q21H",
  # "C64513",  "Q22H",
  # "C64514",  "Q23H",
  # "C64515",  "Q24H",
  # "C64516",  "Q2H",
  "C64536",  "Q2M",
  # "C89791",  "Q36H",
  "C64533",  "Q3D",
  # "C64517",  "Q3H",
  "C64537",  "Q3M",
  # "C139183",  "Q45MIN",
  # "C89790",  "Q48H",
  "C64534",  "Q4D",
  # "C64518",  "Q4H",
  "C64538",  "Q4M",
  "C71124",  "Q5D",
  # "C64519",  "Q5H",
  "C161335",  "Q6D",
  # "C64520",  "Q6H",
  "C89789",  "Q6M",
  # "C174288",  "Q72H",
  "C139177",  "Q7D",
  # "C64521",  "Q7H",
  # "C64523",  "Q8H",
  # "C64524",  "Q9H",
  "C64595",  "QAM",
  "C25473",  "QD",
  # "C64510",  "QH",
  "C64593",  "QHS",
  "C64530",  "QID",
  "C64498",  "QM",
  "C139178",  "QN",
  "C64525",  "QOD",
  "C64525",  "Q2D",
  "C64596",  "QPM",
  "C64527",  "TID"
  ) %>%
mutate(DoseCount = case_when(
    str_detect(CDISCSubmissionValue, "PER") ~
      as.numeric(str_remove_all(CDISCSubmissionValue, "[\\D]")),
    str_detect(CDISCSubmissionValue, "^(Q|EVERY)\\s?\\d{1,2}") ~
      1 / as.numeric(str_remove_all(CDISCSubmissionValue, "[\\D]")),
    str_detect(CDISCSubmissionValue, "^EVERY (A|E|W)[:alpha:]+") ~ 1,
    str_detect(CDISCSubmissionValue, "^Q(AM|PM|M|N|D|HS)|^PA$") ~ 1,
    str_detect(CDISCSubmissionValue, "BI[DM]") ~ 2,
    str_detect(CDISCSubmissionValue, "TID") ~ 3,
    str_detect(CDISCSubmissionValue, "QID") ~ 4,
    str_detect(CDISCSubmissionValue, "QOD") ~ 0.5,
    ),
  DoseWindow = case_when(
    str_detect(CDISCSubmissionValue, "EVERY \\d{1,2}|PER") ~
      str_remove_all(sub(".* (\\w+)$", "\\1", CDISCSubmissionValue), "S"),
    str_detect(CDISCSubmissionValue, "^Q\\d{1,2}D$") ~ "DAY",
    str_detect(CDISCSubmissionValue, "^Q\\d{1,2}M$") ~ "MONTH",
    CDISCSubmissionValue %in% c("EVERY AFTERNOON", "EVERY EVENING") ~ "DAY",
    CDISCSubmissionValue %in% c("EVERY WEEK") ~ "WEEK",
    CDISCSubmissionValue %in% c("BID", "TID", "QAM", "QPM", "QHS",
                                "QD", "QN", "QID", "QOD") ~ "DAY",
    CDISCSubmissionValue %in% c("QM", "BIM") ~ "MONTH",
    CDISCSubmissionValue == "PA" ~ "YEAR",
  )) %>%
mutate(
  DayConversionFactor = case_when(
    DoseWindow == "DAY" ~ 1,
    DoseWindow == "WEEK" ~ (1 / 7),
    DoseWindow == "MONTH" ~ (1 / 30.4375),
    DoseWindow == "YEAR" ~ (1 / 365.25),
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
#' @param dose_freq The dose frequency
#'
#'   The aggregate dosing frequency used for multiple doses in a row.
#'
#'   Default: EXDOSFRQ
#'
#'   Permitted Values: defined by lookup table.
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
#' @param lookup_table The dose frequency value lookup table
#'
#'   The table used to look up `dose_freq` values and determine the appropriate
#'   multiplier to be used for row generation. If a lookup table other than the
#'   default is used, it must have columns `DoseWindow`, `DoseCount`, and
#'   `DayConversionFactor`
#'
#'   Default: `dose_freq_lookup`
#'
#'  Permitted Values for `DoseWindow`: "DAY", "WEEK", "MONTH", "YEAR"
#'
#' @param lookup_column The dose frequency value column in the lookup table
#'
#'   The column of `lookup_table`.
#'
#'   Default: `CDISCSubmissionValue` (column of `dose_freq_lookup`)
#'
#'
#' @details Each aggregate dose row is split into multiple rows which each represent
#' a single dose.The number of completed dose periods between `start_date` and
#' `end_date` is calculated with `compute_duration` and multiplied by the count
#' of doses per period. `start_date` and `end_date` are then adjusted by an
#' appropriate increment using `DayConversionFactor`.
#'
#' @author Michael Thorpe, Andrew Smith
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
                                       end_date = AENDT,
                                       lookup_table = dose_freq_lookup,
                                       lookup_column = CDISCSubmissionValue) {

  col_names <- colnames(dataset)
  by_vars <- assert_vars(by_vars)
  dose_freq <- assert_symbol(enquo(dose_freq))
  lookup_column <- assert_symbol(enquo(lookup_column))
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, required_vars = quo_c(by_vars, dose_freq, start_date, end_date))
  assert_data_frame(lookup_table, required_vars = vars(DoseWindow, DoseCount, DayConversionFactor))

    # Set up lookup table to be joined to dataset

    lookup <- lookup_table %>%
      rename(!!dose_freq := !!lookup_column)

    # Use compute_duration to determine the number of completed dose periods

    dataset <- dataset %>%
      left_join(lookup, by = as.character(quo_get_expr(dose_freq))) %>%
      mutate(dose_periods =
               compute_duration(!!start_date, !!end_date, out_unit = "days") * DayConversionFactor
      ) %>%
      mutate(dose_count = ceiling(dose_periods * DoseCount))

    # Generate a row for each completed dose

    dataset <- dataset[rep(row.names(dataset), dataset$dose_count), ]

    # Determine amount of days to adjust start_date and end_date

    dataset <- dataset %>%
      group_by(!!!by_vars, !!dose_freq, !!start_date, !!end_date) %>%
      mutate(time_increment = (row_number() - 1) / (DoseCount * (DayConversionFactor))) %>%
      ungroup() %>%
      mutate(day_differential = days(floor(time_increment)))

    # Adjust start_date and end_date, drop calculation columns

    dataset <- dataset %>%
      mutate(!!dose_freq := "ONCE",
             !!end_date := !!start_date + day_differential,
             !!start_date := !!start_date + day_differential
      ) %>%
      select(!!!vars(col_names))

  return(dataset)
}
