#' Derive On-Treatment Flag
#'
#' Derive on-treatment flag (`ONTRTFL`) in an ADaM BDS dataset with a single
#' assessment date
#'
#' @param bds_dataset `data.frame`.
#'
#' @param start_date The assessment date (e.g. `ADT` for the date of a VS test)
#'   Required; A date or date-time object column is expected
#'
#' @param ref_start_date The lower bound of the on-treatment period
#'   Required; A date or date-time object column is expected.
#'
#' @param ref_end_date The upper bound of the on-treatment period
#'   A date or date-time object column is expected.
#'   Optional; This can be null everything after `ref_start_date` is considered
#'   on-treatment.
#'
#' @param new_var The variable to be created.
#'   Required; default is `ONTRTFL`
#'   Values: Y or NA
#'
#' @param ref_end_window A window to add to the upper bound `ref_end_date`.
#'   (e.g. 7 if 7 days should be added to the upper bound)
#'   Optional; default is 0
#'   `ref_end_window` and `ref_end_window_units` should be used together
#'
#' @param ref_end_window_units The units of `ref_end_window`
#'   (e.g. 7 days would have `ref_end_window_units` set to days)
#'   Optional; default is "days", valid values: days, weeks, years
#'   `ref_end_window` and `ref_end_window_units` should be used together
#'
#' @param timepoint_var Timepoint variable used to exclude observations
#'   that occur on the `ref_start_date`.
#'   Optional; default is NA
#'   `timepoint_var` and `start_timepoint_pre_value` should be used together
#'
#' @param start_timepoint_pre_value If a baseline timepoint is used to
#'   differentiate observations taken on the day of treatment, those to be
#'   considered not on-treatment should be listed.
#'   Optional; default is NA
#'   Valid values can be specified as "PRE" or c("PRE", "BEFORE")
#'   `timepoint_var` and `start_timepoint_pre_value` should be used together
#'
#'
#' @details
#' On-Treatment is calculated by determining whether the assessment date or
#' start/stop dates fall between 2 dates. The following logic is used to
#' assign on-treatment = Y:
#'   1. `start_date` is missing and `ref_start_date`is non-missing
#'   2. No timepoint is provided (`timepoint_var`) and both `start_date` and
#'      `ref_start_date` are non-missing and `start_date` = `ref_start_date`
#'   3. A timepoint is provided (`timepoint_var`) and both `start_date` and
#'      `ref_start_date` are non-missing and `start_date` = `ref_start_date` and
#'      `timepoint_var` value is not in `start_timepoint_pre_value` list
#'   4. `ref_end_date` is not provided in the function call and `ref_start_date` <
#'      `start_date`
#'   5. `ref_end_date` is provided in the function call and `ref_start_date` <
#'      `start_date` <= `ref_end_date` + `ref_end_window`.
#'
#' Any date imputations needed should be done prior to calling this function.
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with an additional column named by default
#' `ONTRTFL`

#' @export
#'
#' @examples
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL
#'   "P01",    lubridate::ymd("2020-02-24"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P02",    lubridate::ymd("2020-01-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P03",    lubridate::ymd("2019-12-31"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "N")
#' )
#' #' derive_var_ontrtfl(advs, ADT, TRTSDT, TRTEDT)
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL
#'   "P01",    lubridate::ymd("2020-07-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "N",
#'   "P02",    lubridate::ymd("2020-04-30"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P03",    lubridate::ymd("2020-03-15"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y")
#' )
#' derive_var_ontrtfl(advs, ADT, TRTSDT, TRTEDT, ref_end_window=60,
#'   ref_end_window_units="days")
#'
#' advs <- tibble::tribble(
#' ~USUBJID, ~ADTM, ~TRTSDTM, ~TRTEDTM, ~TPT, ~ONTRTFL,
#' "P01",    lubridate::ymd("2020-01-02T12:00"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "", "Y",
#' "P02",    lubridate::ymd("2020-01-01"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "PRE", NA,
#' "P03",    lubridate::ymd("2019-12-31"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "", NA
#' )
#' derive_var_ontrtfl(advs, ADTM, TRTSDTM, TRTEDTM, timepoint_var=TPT,
#'   start_timepoint_pre_value="PRE")
#'
derive_var_ontrtfl <- function(bds_dataset,
                               start_date,
                               ref_start_date,
                               ref_end_date,
                               new_var = ONTRTFL,
                               ref_end_window = 0,
                               ref_end_window_units = "days",
                               timepoint_var = NA,
                               start_timepoint_pre_value = NA) {

   new_var <- enquo(new_var)
   start_date <- enquo(start_date)
   ref_start_date <- enquo(ref_start_date)
   ref_end_date <- enquo(ref_end_date)
   ref_end_window <- enquo(ref_end_window)
   timepoint_var <- enquo(timepoint_var)
   start_timepoint_pre_value = enquo(start_timepoint_pre_value)

   '%!in%' <- Negate(`%in%`)

   #Check validity of window offset and units
   if (ref_end_window_units %in% c("days", "weeks", "years")) {
      if (is.na(as.integer(quo_text(ref_end_window))) == TRUE) {
         msg <- paste0("The window is missing or invalid [ref_end_window]. Set to valid number or remove units if not applicable")
         abort(msg)
      }
   } else if(ref_end_window_units %!in% c("days", "weeks", "years") &
             as.integer(quo_text(ref_end_window)) != 0) {
      msg <- paste0("The window unit [ref_end_window_units] is not days, weeks, years and the interval [ref_end_window] is specified")
      abort(msg)
   }

   #Check to ensure dates are passed

   #Check if timepoint passed, it is a valid column in the dataframe

   #Create a temporary ref_end_date for computations that includes the
   #ref_end_window
   if (quo_text(ref_end_date) == "") {
      ontrtfl <- bds_dataset %>%
         mutate(TEMP_REF_END_DATE = NA)
   }else {
      ontrtfl <- bds_dataset %>%
         mutate(TEMP_REF_END_DATE =
                   case_when(!!ref_end_window_units == "days" ~
                                !!ref_end_date + days(x = !!ref_end_window),
                             !!ref_end_window_units == "weeks" ~
                                !!ref_end_date + weeks(x = !!ref_end_window),
                             !!ref_end_window_units == "years" ~
                                !!ref_end_date + years(x = !!ref_end_window),
                             is.na(!!ref_end_date) == FALSE ~ !!ref_end_date))
   }

   #Derive On-Treatment flag
   #Scenario 1: No treatment end date is passed
   #Scenario 2: Treatment end date is passed, window added above
   ontrtfl <- ontrtfl %>%
      mutate(!!new_var :=
                case_when(is.na(!!start_date) == TRUE &
                             is.na(!!ref_start_date) == FALSE ~ "Y",
                          is.na(!!timepoint_var) == TRUE &
                             is.na(!!start_date) == FALSE &
                             is.na(!!ref_start_date) == FALSE &
                             !!ref_start_date == !!start_date ~ "Y",
                          is.na(!!timepoint_var) == FALSE &
                             is.na(!!start_date) == FALSE &
                             is.na(!!ref_start_date) == FALSE &
                             !!ref_start_date == !!start_date &
                             !!timepoint_var %!in% !!start_timepoint_pre_value
                          ~ "Y"))

   if (quo_text(ref_end_date) == "") {
      ontrtfl <- ontrtfl %>%
         mutate(!!new_var :=
                   ifelse(is.na(!!ref_start_date) == FALSE &
                             !!ref_start_date < !!start_date, "Y", !!new_var))
   } else {
      ontrtfl <- ontrtfl %>%
         mutate(!!new_var :=
                   ifelse(is.na(!!ref_start_date) == FALSE &
                             !!ref_start_date < !!start_date &
                             is.na(TEMP_REF_END_DATE) == FALSE &
                             !!start_date <= TEMP_REF_END_DATE, "Y", !!new_var))
   }

   ontrtfl <- ontrtfl %>%
      select(-TEMP_REF_END_DATE)

   ontrtfl
}
