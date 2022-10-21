admiral_options <- list(subject_keys = vars(STUDYID, USUBJID),
                        future_input = vars(STUDYID, SITEID))
#' Call a package-standard input
#'
#' Call a package-standard input that can be modified for advanced users.
#'
#' @param input The function input necessary.
#'
#' @author Zelos Zhu
#'
#' @return
#' Call on an admiral option to be called for function inputs: e.g `subject_keys`.
#'
#' @family XXX
#' @keywords XXX
#'
#' @export
#'
#' @seealso [vars()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' get_admiral_options(subject_keys)
get_admiral_options <- function(input){
  #Check for valid input - catch function abuse
  assert_expr(enquo(input))

  #Find which admiral_options is being called upon
  x = enquo(input)
  if(as_name(x) %in% names(admiral_options)) {
      index = which(as_name(x) == names(admiral_options))
      return(admiral_options[[index]])
  }

  #Return message otherwise, catch typos
  else {
    default_err_msg <- sprintf(paste("Invalid function argument, select one unquoted of:",
                                     paste0(names(admiral_options), collapse = " or ")
                                     )
                               )
    abort(default_err_msg)
  }
}

#' Set a package-standard input
#'
#' Set a package-standard input that can be modified for advanced users.
#'
#' @param subject_keys Subjects used for several `derive_()` functions, defaults to
#'   vars(STUDYID, USUBJID)
#'
#' @param future_input Possible future input to figure out how to scale this
#'
#' @author Zelos Zhu
#'
#' @return
#' Modify an admiral option which will affect downstream function inputs: e.g `subject_keys`.
#'
#' @family XXX
#' @keywords XXX
#'
#' @export
#'
#' @seealso [vars()]
#'
#' @examples
#' set_admiral_options(subject_keys = vars(STUDYID, USUBJID2))
set_admiral_options <- function(subject_keys,
                                future_input){
  if(!missing(subject_keys)){
    assert_vars(subject_keys)
    admiral_options$subject_keys <<- subject_keys
  }

  if(!missing(future_input)){
    assert_vars(future_input)
    admiral_options$future_input <<- future_input
  }
}

