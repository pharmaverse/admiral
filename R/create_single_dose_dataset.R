#' Pre-Defined Dose Frequencies
#'
#' These pre-defined dose frequencies are sourced from
#' [CDISC](https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.pdf). The
#' number of rows to generate using `create_single_dose_dataset()` arguments
#' `start_date` and `end_date` is derived from `DOSE_COUNT`, `DOSE_WINDOW`, and
#' `CONVERSION_FACTOR` with appropriate functions from `lubridate`.
#'
#' @details
#' `NCI_CODE` and `CDISC_VALUE` are included from the CDISC source for
#' traceability.
#'
#' `DOSE_COUNT` represents the number of doses received in one single unit of
#' `DOSE_WINDOW`. For example, for `CDISC_VALUE=="10 DAYS PER MONTH`,
#' `DOSE_WINDOW=="MONTH"` and `DOSE_COUNT==10`. Similarly, for
#' `CDISC_VALUE=="EVERY 2 WEEKS`, `DOSE_WINDOW=="WEEK"` and
#' `DOSE_COUNT==0.5` (to yield one dose every two weeks).
#'
#' `CONVERSION_FACTOR` is used to convert `DOSE_WINDOW` units `"WEEK"`,
#'  `"MONTH"`, and `"YEAR"` to the unit `"DAY"`.
#'
#' For example, for `CDISC_VALUE=="10 DAYS PER MONTH`, `CONVERSION_FACTOR`
#' is `0.0329`. One day of a month is assumed to be `1 / 30.4375` of a month (one
#' day is assumed to be `1/365.25` of a year).
#' Given only `start_date` and `end_date` in the aggregate dataset, `CONVERSION_FACTOR`
#' is used to calculate specific dates for`start_date` and `end_date` in the
#' resulting single dose dataset for the doses that occur. In such cases, doses
#' are assumed to occur at evenly spaced increments over the interval.
#'
#'
#' To see the entire table in the console, run `print(dose_freq_lookup)`.
#'
#' @seealso [create_single_dose_dataset()]
#'
#' @export
#'
#' @keywords metadata
#'
#' @rdname dose_freq_lookup
#'
dose_freq_lookup <- tibble::tribble(
  ~NCI_CODE, ~CDISC_VALUE,
  "C64526", "1 TIME PER WEEK",
  "C139179", "10 DAYS PER MONTH",
  "C64497", "2 TIMES PER WEEK",
  "C98861", "2 TIMES PER YEAR",
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
  "C64500",  "Q10H",
  "C64501",  "Q11H",
  "C64502",  "Q12H",
  "C64503",  "Q13H",
  "C64504",  "Q14H",
  "C64505",  "Q15H",
  "C64506",  "Q16H",
  "C64507",  "Q17H",
  "C64508",  "Q18H",
  "C64509",  "Q19H",
  "C64511",  "Q20H",
  "C64512",  "Q21H",
  "C64513",  "Q22H",
  "C64514",  "Q23H",
  "C64515",  "Q24H",
  "C64516",  "Q2H",
  "C64536",  "Q2M",
  "C89791",  "Q36H",
  "C64533",  "Q3D",
  "C64517",  "Q3H",
  "C64537",  "Q3M",
  "C139183",  "Q45MIN",
  "C89790",  "Q48H",
  "C64534",  "Q4D",
  "C64518",  "Q4H",
  "C64538",  "Q4M",
  "C71124",  "Q5D",
  "C64519",  "Q5H",
  "C161335",  "Q6D",
  "C64520",  "Q6H",
  "C89789",  "Q6M",
  "C174288",  "Q72H",
  "C139177",  "Q7D",
  "C64521",  "Q7H",
  "C64523",  "Q8H",
  "C64524",  "Q9H",
  "C64595",  "QAM",
  "C25473",  "QD",
  "C64510",  "QH",
  "C64593",  "QHS",
  "C64530",  "QID",
  "C64498",  "QM",
  "C139178",  "QN",
  "C64525",  "QOD",
  "C64525",  "Q2D",
  "C64596",  "QPM",
  "C64527",  "TID"
  ) %>%
mutate(DOSE_COUNT = case_when(
    str_detect(CDISC_VALUE, "PER [WMY]") ~
      as.numeric(str_remove_all(CDISC_VALUE, "[\\D]")),
    str_detect(CDISC_VALUE, "PER [D]") ~
      24 / as.numeric(str_remove_all(CDISC_VALUE, "[\\D]")),
    str_detect(CDISC_VALUE, "^Q\\d{1,2}(H|MIN)") ~
      1 / as.numeric(str_remove_all(CDISC_VALUE, "[\\D]")),
    str_detect(CDISC_VALUE, "^(Q|EVERY)\\s?\\d{1,2}") ~
      1 / as.numeric(str_remove_all(CDISC_VALUE, "[\\D]")),
    str_detect(CDISC_VALUE, "^EVERY (A|E|W)[:alpha:]+") ~ 1,
    str_detect(CDISC_VALUE, "^Q(AM|PM|M|N|D|HS)|^PA$") ~ 1,
    str_detect(CDISC_VALUE, "^QH$") ~ 1,
    str_detect(CDISC_VALUE, "BI[DM]") ~ 2,
    str_detect(CDISC_VALUE, "TID") ~ 3,
    str_detect(CDISC_VALUE, "QID") ~ 4,
    str_detect(CDISC_VALUE, "QOD") ~ 0.5,
    ),
  DOSE_WINDOW = case_when(
    str_detect(CDISC_VALUE, "EVERY \\d{1,2}|PER [WMY]") ~
      str_remove_all(sub(".* (\\w+)$", "\\1", CDISC_VALUE), "S"),
    str_detect(CDISC_VALUE, "^Q\\d{1,2}D$") ~ "DAY",
    str_detect(CDISC_VALUE, "^Q\\d{1,2}M$") ~ "MONTH",
    str_detect(CDISC_VALUE, "^Q\\d{0,2}H$|PER D") ~ "HOUR",
    str_detect(CDISC_VALUE, "^Q\\d{1,2}MIN$") ~ "MINUTE",
    CDISC_VALUE %in% c("EVERY AFTERNOON", "EVERY EVENING") ~ "DAY",
    CDISC_VALUE %in% c("EVERY WEEK") ~ "WEEK",
    CDISC_VALUE %in% c("BID", "TID", "QAM", "QPM", "QHS",
                                "QD", "QN", "QID", "QOD") ~ "DAY",
    CDISC_VALUE %in% c("QM", "BIM") ~ "MONTH",
    CDISC_VALUE == "PA" ~ "YEAR",
  )) %>%
mutate(
  CONVERSION_FACTOR = case_when(
    DOSE_WINDOW == "MINUTE" ~ 1,
    DOSE_WINDOW == "HOUR" ~ 1,
    DOSE_WINDOW == "DAY" ~ 1,
    DOSE_WINDOW == "WEEK" ~ (1 / 7),
    DOSE_WINDOW == "MONTH" ~ (1 / 30.4375),
    DOSE_WINDOW == "YEAR" ~ (1 / 365.25),
  )
)

#' Create dataset of single doses
#'
#' Derives dataset of single dose from aggregate dose information. This may be
#' necessary when e.g. calculating last dose before an adverse event in `ADAE`
#' or deriving a total dose parameter in `ADEX` when `EXDOSFRQ != ONCE`.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by `dose_freq`, `start_date` and the `end_date`
#'   parameters are expected.
#'
#' @param dose_freq The dose frequency
#'
#'   The aggregate dosing frequency used for multiple doses in a row.
#'
#'   Default: `EXDOSFRQ`
#'
#'   Permitted Values: defined by lookup table.
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_var_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Default: `ASTDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_var_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Default: `AENDT`
#'
#' @param lookup_table The dose frequency value lookup table
#'
#'   The table used to look up `dose_freq` values and determine the appropriate
#'   multiplier to be used for row generation. If a lookup table other than the
#'   default is used, it must have columns `DOSE_WINDOW`, `DOSE_COUNT`, and
#'   `CONVERSION_FACTOR`. The default table `dose_freq_lookup` is described in
#'   detail [here][dose_freq_lookup].
#'
#'   Default: `dose_freq_lookup`
#'
#'   Permitted Values for `DOSE_WINDOW`: `"MINUTE"`, `"HOUR"`, `"DAY"`,
#'   `"WEEK"`, `"MONTH"`, `"YEAR"`
#'
#' @param lookup_column The dose frequency value column in the lookup table
#'
#'   The column of `lookup_table`.
#'
#'   Default: `CDISC_VALUE` (column of `dose_freq_lookup`)
#'
#'
#' @details Each aggregate dose row is split into multiple rows which each
#'   represent a single dose.The number of completed dose periods between
#'   `start_date` and `end_date` is calculated with `compute_duration` and
#'   multiplied by `DOSE_COUNT`. For `DOSE_WINDOW` values of `"WEEK"`,
#'   `"MONTH"`, and `"YEAR"`, `CONVERSION_FACTOR` is used to convert into days
#'   the time object to be added to `start_date`.
#'
#' @author Michael Thorpe, Andrew Smith
#'
#' @keywords adae adex user_utility
#'
#' @return The input dataset with a single dose per row.
#'
#' @export
#'
#' @examples
#' # Example with default lookup
#'
#' library(lubridate)
#'
#' data <- tibble::tribble(
#'   ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
#'   "P01", "Q2D", ymd("2021-01-01"), ymd("2021-01-07"),
#'   "P01", "Q3D", ymd("2021-01-08"), ymd("2021-01-15"),
#'   "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd("2021-01-29")
#' )
#'
#' create_single_dose_dataset(data)
#'
#' # Example with custom lookup
#'
#' custom_lookup <- tibble::tribble(
#'   ~Value,  ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
#'   "Q30MIN", (1 / 30), "MINUTE", 1,
#'   "Q90MIN", (1 / 90), "MINUTE", 1
#' )
#'
#' data <- tibble::tribble(
#'   ~USUBJID, ~EXDOSFRQ, ~ASTDTM, ~AENDTM,
#'   "P01", "Q30MIN", ymd_hms("2021-01-01T06:00:00"), ymd_hms("2021-01-01T07:00:00"),
#'   "P02", "Q90MIN", ymd_hms("2021-01-01T06:00:00"), ymd_hms("2021-01-01T09:00:00")
#' )
#'
#' create_single_dose_dataset(data,
#'                            lookup_table = custom_lookup,
#'                            lookup_column = Value,
#'                            start_date = ASTDTM,
#'                            end_date = AENDTM)

create_single_dose_dataset <- function(dataset,
                                       dose_freq = EXDOSFRQ,
                                       start_date = ASTDT,
                                       end_date = AENDT,
                                       lookup_table = dose_freq_lookup,
                                       lookup_column = CDISC_VALUE) {

  col_names <- colnames(dataset)
  dose_freq <- assert_symbol(enquo(dose_freq))
  lookup_column <- assert_symbol(enquo(lookup_column))
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, required_vars = quo_c(dose_freq, start_date, end_date))
  assert_data_frame(lookup_table, required_vars = vars(DOSE_WINDOW, DOSE_COUNT, CONVERSION_FACTOR))

  # Set up lookup table to be joined to dataset

  lookup <- lookup_table %>%
    rename(!!dose_freq := !!lookup_column)

  # Check values of lookup vs. data and return error if values are not covered

  value_check <- dataset %>%
    select(!!dose_freq) %>%
    anti_join(lookup, by = as.character(quo_get_expr(dose_freq))) %>%
    unique()

  if (nrow(value_check) > 0) {
    values_not_found <- paste0(value_check %>% select(!!dose_freq), collapse = ", ")

    err_msg <- paste0(
      sprintf("The following values of %s in %s do not appear in %s:\n",
              as.character(quo_get_expr(dose_freq)),
              arg_name(substitute(dataset)),
              arg_name(substitute(lookup_table))
              ), values_not_found)
    abort(err_msg)
  }

  # Use compute_duration to determine the number of completed dose periods

  dataset_part_1 <- dataset %>%
    filter(!!dose_freq == "ONCE")

  dataset_part_2 <- dataset %>%
    filter(!!dose_freq != "ONCE")

  dataset_part_2 <- dataset_part_2 %>%
    left_join(lookup, by = as.character(quo_get_expr(dose_freq))) %>%
    mutate(dose_periods = case_when(
        DOSE_WINDOW == "MINUTE" ~ compute_duration(!!start_date, !!end_date, out_unit = "minutes"),
        DOSE_WINDOW == "HOUR" ~ compute_duration(!!start_date, !!end_date, out_unit = "hours"),
        DOSE_WINDOW == "DAY" ~ compute_duration(!!start_date, !!end_date, out_unit = "days"),
        DOSE_WINDOW == "WEEK" ~ compute_duration(!!start_date, !!end_date, out_unit = "weeks"),
        DOSE_WINDOW == "MONTH" ~ compute_duration(!!start_date, !!end_date, out_unit = "months"),
        DOSE_WINDOW == "YEAR" ~ compute_duration(!!start_date, !!end_date, out_unit = "years")
    )
      ) %>%
    mutate(dose_count = ceiling(dose_periods * DOSE_COUNT)) %>%
    derive_var_obs_number(new_var = grpseq)

  # Flag to determine if date or datetime must be returned

  time_flag <- nrow(dataset_part_2 %>%
                       filter(DOSE_WINDOW %in% c("MINUTE", "HOUR"))) > 0


  if (time_flag &
      (is.Date(eval_tidy(start_date, dataset)) | is.Date(eval_tidy(end_date, dataset)))) {
    err_msg <- paste0(
      sprintf("%s involves hours or minutes but one of %s or %s is a date variable. ",
              as.character(quo_get_expr(dose_freq)),
              as.character(quo_get_expr(start_date)),
              as.character(quo_get_expr(end_date))),
      "\nPlease provide datetime variables for start_date and end_date arguments.")
    abort(err_msg)
  }

  # Generate a row for each completed dose

  dataset_part_2 <- dataset_part_2[rep(row.names(dataset_part_2), dataset_part_2$dose_count), ]

  # Determine amount of days to adjust start_date and end_date

  dataset_part_2 <- dataset_part_2 %>%
    group_by(grpseq, !!dose_freq, !!start_date, !!end_date) %>%
    mutate(time_increment = (row_number() - 1) / (DOSE_COUNT)) %>%
    ungroup() %>%
    mutate(time_differential = case_when(
      DOSE_WINDOW == "MINUTE" ~ minutes(floor(time_increment)),
      DOSE_WINDOW == "HOUR" ~ hours(floor(time_increment)),
      DOSE_WINDOW %in% c("DAY", "WEEK", "MONTH", "YEAR") ~
        days(floor(time_increment / CONVERSION_FACTOR))
           )
    )

  # Adjust start_date and end_date, drop calculation columns, make sure nothing
  # later than end_date shows up in output

  dataset_part_2 <- dataset_part_2 %>%
    mutate(!!dose_freq := "ONCE",
    !!start_date := !!start_date + time_differential)

  dataset_part_2 <- dataset_part_2 %>%
    filter(!(!!start_date > !!end_date)) %>%
    mutate(!!end_date := !!start_date) %>%
    select(!!!vars(all_of(col_names)))

  #Stitch back together

  dataset <- bind_rows(dataset_part_1, dataset_part_2)

  return(dataset)
}
