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
      " argument",
      if_else(length(expected_vars) > 1, "s", ""),
      " are expected to be in the dataset."
    )
  }
  return(dataset_text)
}


# function to properly document the by_vars argument including specifying the use of exprs():

roxygen_param_by_vars <- function(rename = FALSE) {
  by_vars_text <- ""

  if (rename) {
    by_vars_text <- paste0(
      by_vars_text,
      "Variables can be renamed by naming the element, i.e. \n",
      "`by_vars = exprs(<name in input dataset> = <name in additional dataset>)`, ",
      "similar to the `dplyr` joins.\n \n"
    )
  }

  by_vars_text <- paste0(
    by_vars_text,
    "*Permitted Values*: list of variables created by `exprs()` \n",
    "e.g. `exprs(USUBJID, VISIT)`"
  )

  return(by_vars_text)
}

roxygen_order_na_handling <- function() {
  paste(
    "For handling of `NA`s in sorting variables see",
    "[Sort Order](../articles/generic.html#sort_order)."
  )
}
