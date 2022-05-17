#' Calls a Function Provided by the User
#'
#' Calls a function provided by the user and adds the function call to the error
#' message if the call fails.
#'
#' @param call Call to be executed
#'
#' @author Stefan Bundfuss
#'
#' @return The return value of the function call
#'
#' @keywords dev_utility
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
  tryCatch(
    eval_tidy(call),
    error = function(cnd) {
      abort(
        paste0("Calling ", rlang::as_label(enexpr(call)), " caused the following error:\n", cnd)
      )
    }
  )
}
