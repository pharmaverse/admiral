# THIS FILE INCLUDES FUNCTION THAT ARE RE-USED THROUGHOUT
# roxygen2 COMMENTS TO GENERATE DOCUMENTATION TEXT

roxygen_param_dataset <- function(expected_vars = NULL) {
  if (is.null(expected_vars)) {
    dataset_text <- "Input dataset"
  } else {
    dataset_text <- paste0(
      "Input dataset \n \n",
      "The variables specified by the ",
      enumerate(expected_vars),
      " argument(s) to be expected."
    )
  }
  return(dataset_text)
}


# function to properly document the by_vars argument including specifying the use of exprs():

roxygen_param_by_vars <- function(additional_dataset = NULL, unique = FALSE, rename = FALSE) {
  by_vars_text <- ""

  if (!is.null(additional_dataset)) {
    if (unique) {
      by_vars_text <- paste0(
        by_vars_text,
        "Variables must be a unique key of the selected observations in ",
        enumerate(additional_dataset), ". \n \n"
      )
    }
    if (rename) {
      by_vars_text <- paste0(
        by_vars_text,
        "Variables from ", enumerate(additional_dataset),
        " can be renamed by naming the element, I.e. \n",
        "by_vars = exprs(<name in input dataset> = <name in additional dataset>),",
        "similar to the dplyr joins.\n \n"
      )
    }
  }

  by_vars_text <- paste0(
    by_vars_text,
    "*Permitted Values*: list of variables created by exprs() \n \n",
    "e.g. exprs(USUBJID, VISIT)"
  )

  return(by_vars_text)
}
