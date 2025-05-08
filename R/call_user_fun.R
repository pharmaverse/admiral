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
    details = ""
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
