call_user_fun <- function(call) {
  tryCatch(
    call,
    error = function(cnd) {
      abort(
          paste0("Calling ", rlang::as_label(enexpr(call)), " caused the following error:\n", cnd)
      )
    }
  )
}
