#' Calls a Function Provided by the User
#'
#' @description
#' r lifecycle::badge("deprecated")
#' Calls a function provided by the user and adds the function call to the error
#' message if the call fails.

#' @param call Call to be executed
#'
#'
#' @return The return value of the function call
#'
#' @family  deprecated
#' @keywords deprecated
#'
#' @export
#'
#' @examples
#' call_user_fun(compute_bmi(
#'   height = 172,
#'   weight = 60
#' ))
#'
#' try(call_user_fun(compute_bmi(
#'   height = 172,
#'   weight = "hallo"
#' )))
call_user_fun <- function(call) {
  deprecate_inform(
    when = "1.3.0",
    what = "call_user_fun()",
    details = c(
      "`call_user_fun()` is no longer supported and no replacement is provided; ",
      "The original code for this function is here: ",
      "https://github.com/pharmaverse/admiral/blob/v1.2.0/R/call_user_fun.R#L26-L39"
    )
  )
  tryCatch(
    eval_tidy(call),
    error = function(cnd) {
      cli_abort(
        message = c(
          "Calling {.code {as_label(enexpr(call))}} caused the following error:",
          conditionMessage(cnd)
        ),
        call = parent.frame(n = 4)
      )
    }
  )
}
