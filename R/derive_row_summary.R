#' Add new row within by groups using aggregation functions
#'
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple rows. The ADaM basic dataset structure
#' variable `DTYPE` is available to indicate when a new derived row has been added to a
#' dataset.
#'
#' @param dataset A data frame.
#' @param by_vars Variables to consider for generation of groupwise summary
#'   rows. Providing the names of variables in [c()] will create a groupwise
#'   summary and generate summary rows for the specified groups.
#' @param fns List of formulas specifying variable to use for aggregations.
#'   This can include base functions like `mean`, `min`, `max`, `median`, `sd`,
#'   or `sum` or any other user-defined aggregation function. For example,
#'   `fns = list(AVAL ~ mean, CHG ~ sum(., na.rm = TRUE)`. In the formula
#'   representation, a `.` serves as the data to be summarized
#'   (e.g., `sum(CHG, na.rm = TRUE)`).
#' @param filter Logical expression indicating rows to keep during aggregations.
#' @param values_set A list of variable name-value pairs. Use this argument if
#'   you need to change the values of any newly derived rows. For example,
#'   `values_set = list(AVISITN = 9999, AVISIT= "Endpoint")` would change of
#'   value of AVISITN to 9999 and AVISIT to Endpoint instead of retaining.
#' @param values_drop Providing the names of variables in [c()] will drop values
#'   and set as missing.
#'
#' @return
#' @export
#'
#' @examples
#' # In this example, the analysis requirement is to summarize the average of
#' # the triplicate ECG interval values (AVAL).
#'
#' adeg <- tibble::tribble(
#' ~USUBJID,  ~EGSEQ,~EGREPNUM,~PARAM,~VISIT,~AVISIT,~EGDTC,~AVAL,~TRTA,~SAFFL,
#' "XYZ-1001",1,1,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:50:16",385,"","Y",
#' "XYZ-1001",2,2,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:52:59",399,"","Y",
#' "XYZ-1001",3,3,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:56:07",396,"","Y",
#' "XYZ-1001",4,1,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:45:11",384,"Placebo","Y",
#' "XYZ-1001",5,2,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:48:07",393,"Placebo","Y",
#' "XYZ-1001",6,3,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:51:04",388,"Placebo","Y",
#' "XYZ-1001",7,1,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:45:03",385,"Placebo","Y",
#' "XYZ-1001",8,2,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:48:07",394,"Placebo","Y",
#' "XYZ-1001",9,3,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:51:05",402,"Placebo","Y",
#' "XYZ-1002",1,1,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T07:58:05",399,"","Y",
#' "XYZ-1002",2,2,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T07:58:05",410,"","Y",
#' "XYZ-1002",3,3,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T08:01:06",392,"","Y",
#' "XYZ-1002",4,1,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:50:04",401,"Active 20mg","Y",
#' "XYZ-1002",5,2,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:53:51",407,"Active 20mg","Y",
#' "XYZ-1002",6,3,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:56:21",400,"Active 20mg","Y",
#' "XYZ-1002",7,1,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:50:07",412,"Active 20mg","Y",
#' "XYZ-1002",8,2,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:53:08",414,"Active 20mg","Y",
#' "XYZ-1002",9,3,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:56:05",402,"Active 20mg","Y",
#' )
#' derive_row_summary(adeg,
#'                    by_vars = c(USUBJID, PARAM, VISIT, AVISIT),
#'                    fns = list(AVAL ~ mean(., na.rm = TRUE)),
#'                    values_drop = c(EGSEQ, EGREPNUM, EGDTC))
derive_row_summary <- function(dataset,
                               by_vars,
                               fns,
                               filter = NULL,
                               values_set = list(),
                               values_drop = NULL
                               ) {
  fns <- fns
  group_vars <- dplyr::select(dataset, !! rlang::enquo(by_vars)) %>% names()
  drop_vars <- dplyr::select(dataset, !! rlang::enquo(values_drop)) %>% names()
  browser()

}

# Helper ------------------------------------------------------------------

setup_fns <- function(fns, env = rlang::caller_env()) {
  if (is.function(fns)) {
    list(fun = fns, var = "AVAL")
  } else if (rlang::is_formula(fns)) {
    var <- rlang::f_rhs(fns)[[2]] %>% as_string()
    list(fun = rlang::as_closure(fns), var = var)
  } else if (is.list(fns)) {
    chk_fns <- vapply(fns,
                          function(x) { rlang::is_formula(x, lhs = FALSE) },
                          logical(1L))
    if (length(fns[!chk_fns]) > 0) {
      rlang::abort(
        c(paste("Problem with input `fns` at index", enumerate(which(chk_fns))),
        i = "`fns` must be list of formulas.")
      )
    }

    funs <- lapply(fns, rlang::as_closure)
    vars <- lapply(fns, function(x) {
      rlang::f_rhs(x)[[2]] %>% as_string()
    })
    list(fun = funs, var = vars)
  } else {
    rlang::abort(
      c("Problem with input `fns`.",
      i = "`.fns` must be a function, a formula, or a list of functions/formulas."
    ))
  }
}



