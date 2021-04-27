#' Derive Analysis Baseline Flag
#'
#' Adds the analysis baseline flag (`ABLFL`) to the dataset.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param order Sort order
#'
#'   The first or last observation is determined with respect to the specified
#'   order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param mode Flag mode
#'
#'   Determines of the first or last observation is flagged.
#'
#'   Default: `"last"`
#'
#'   Permitted Values: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param flag_filter Filter for flag data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the parameter is not specified, all observations are
#'   considered.
#'
#'   Permitted Values: a condition
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   Default: `"warning"`
#'
#'   Permitted Values: `"none"`, `"warning"`, `"error"`
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the flag mode) will be flagged
#'   as `ABLFL` = Y in the output dataset.
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' input <- tibble::tribble(
#' ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~DTYPE,
#' "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,
#' "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
#' "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
#' "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
#' "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,
#'
#' "TEST01", "PAT02",  "PARAM01", "SCREENING",as.Date("2021-04-27"), 15.0, "AVERAGE",
#' "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE",
#' "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
#' "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
#' "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE",
#'
#' "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, "AVERAGE",
#' "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-25"), 14.0, "AVERAGE",
#' "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-23"), 15.0, NA,
#' "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "AVERAGE",
#' "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,
#'
#' "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, NA,
#' "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
#' "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
#' "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
#' "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA
#' )
#'
#' Last observation
#' derive_var_ablfl(input,
#'                  by_vars = exprs(USUBJID, PARAMCD),
#'                  order = exprs(ADT),
#'                  flag_filter = expr(AVISIT == "BASELINE"))
#'
#' Worst observation - Direction = High
#' derive_var_ablfl(input,
#'                  by_vars = exprs(USUBJID, PARAMCD),
#'                  order = exprs(AVAL, ADT),
#'                  flag_filter = expr(AVISIT == "BASELINE"))
#'
#' Worst observation - Direction = Lo
#' derive_var_ablfl(input,
#'                  by_vars = exprs(USUBJID, PARAMCD),
#'                  order = exprs(desc(AVAL), ADT),
#'                  flag_filter = expr(AVISIT == "BASELINE"))
#'
#' Average observation
#' derive_var_ablfl(input,
#'                  by_vars = exprs(USUBJID, PARAMCD),
#'                  order = exprs(ADT, desc(AVAL)),
#'                  flag_filter = expr(AVISIT == "BASELINE" &
#'                    DTYPE == "AVERAGE" & is.na(DTYPE) == FALSE))
#'
#' data("vs")

derive_var_ablfl <- function(dataset,
                             by_vars,
                             order,
                             mode = "last",
                             flag_filter,
                             check_type = "warning") {

  dataset <- derive_extreme_flag(dataset,
                                 new_var = ABLFL,
                                 by_vars,
                                 order,
                                 mode = mode,
                                 flag_filter,
                                 check_type = check_type)

}
