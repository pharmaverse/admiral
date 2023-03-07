# Sometimes functions from other libraries become critically important while using admiral
# and thus should be included as an export. This applies especially to functions which are
# frequently called within `admiral` function arguments. The goal of these exports is such that
# admiral comes ready "out of the box", similar to how one might think the pipe operator, `%>%`,
# comes from `dplyr` but is actually native to `magrittr`.

#' rlang exprs
#'
#' See \code{rlang::\link[rlang:exprs]{exprs}} for details.
#'
#' @name exprs
#' @rdname reexport-exprs
#' @keywords reexport
#' @importFrom rlang exprs
#' @export
NULL

#' Create List of Quosures
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `exprs()` instead.
#'
#' @param ... List of variables
#'
#' @return List of expressions
#'
#' @keywords deprecated
#' @family deprecated
#'
#' @export
vars <- function(...) {
  deprecate_warn(
    "0.10.0",
    "vars()",
    "exprs()",
    details = paste(
      "The admiral functions no longer expects lists of quosures created by `vars()`",
      "but lists of expressions created by `exprs()`.",
      "Please update your function calls ASAP."
    )
  )
  exprs(...)
}

# Force admiral definition of vars().
# This and the definition of vars() above should be removed in the release after admiral 0.10.0
setHook(
  packageEvent("dplyr", "attach"), function(...) {
    if (get_admiral_option("force_admiral_vars")) {
      inform(paste0(
        "admiral definition of `vars()` is forced.\n",
        "If you want to use the dplyr definition of `vars()`, call ",
        "`set_admiral_options(force_admiral_vars = FALSE)` ",
        "before attaching dplyr."
      ))
      detach("package:admiral")
      suppressWarnings(library("admiral",
        pos = 2L,
        warn.conflicts = FALSE,
        quietly = TRUE,
        character.only = TRUE,
        verbose = FALSE
      ))
    }
  },
  "append"
)

#' dplyr desc
#'
#' See \code{dplyr::\link[dplyr:desc]{desc}} for details.
#'
#' @name desc
#' @rdname reexport-desc
#' @keywords reexport
#' @importFrom dplyr desc
#' @export
NULL

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:%>%]{%>%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @keywords reexport
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
