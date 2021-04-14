#' Derive On-Treatment Flag
#'
#' Derive on-treatment flag (`ONTRTFL`) in an ADaM dataset
#'
#' @param bds_dataset `data.frame`.
#'   The column specified by the `start_date` is required, default: `ADT`
#'   The column specified by the `ref_start_date'`is required, default: `TRTSDT`
#'   The column specified by the `ref_end_date` is required, default: `TRTEDT`
#'
#' @details
#' On-Treatment is calculated by determining whether the assessment date or
#' start/stop dates fall between 2 dates.
#' If the `start_date` is missing then `ONTRTFL` will be set to null.
#' If the `ref_start_date` is missing and the `ref_end_date` is present and
#' `start_date` is <= `ref_end_date`, then the `ONTRTFL` = Y.
#' If the `ref_end_date` is missing and the `ref_start_date` is present and
#' `start_date` is >= `ref_start_date`, then the `ONTRTFL` = Y.
#' If the both the `ref_start_date` and `ref_end_date` are missing then
#' `ONTRTFL` will be null.
#' Any date imputations needed should be done prior to calling this function.
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with an additional column named `ONTRTFL`
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
#'             lubridate::ymd("2020-03-01"), "N",
#'   "P04",    lubridate::ymd("2020-03-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P05",    lubridate::ymd("2020-04-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "N"
#'   "P06",    lubridate::ymd("2020-01-01"), lubridate::ymd(""),
#'             lubridate::ymd("2020-03-01"), "Y"
#'   "P07",    lubridate::ymd("2020-04-01"), lubridate::ymd(""),
#'             lubridate::ymd("2020-03-01"), "N"
#'   "P08",    lubridate::ymd("2020-03-01"), lubridate::ymd(""),
#'             lubridate::ymd("2020-03-01"), "N"
#'   "P09",    lubridate::ymd("2019-12-31"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd(""), "N"
#'   "P10",    lubridate::ymd("2019-01-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd(""), "Y"
#'   "P11",    lubridate::ymd("2019-01-010"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd(""), "Y"
#'   "P12",    lubridate::ymd(""), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd(""), NA
#'   "P13",    lubridate::ymd("2020-01-01"), lubridate::ymd(""),
#'             lubridate::ymd(""), NA
#' )
#' derive_var_ontrtfl(advs)
#'

derive_var_ontrtfl <- function(bds_dataset, start_date = ADT,
                               ref_start_date = TRTSDT, ref_end_date = TRTEDT) {

}
