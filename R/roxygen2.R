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

roxygen_param_by_vars <- function(type = "group", additional_dataset = NULL) {

  if (type == "group") {
    by_vars_text <- "A list of quoted expressions specifying variables to group by. \n \n"
  }


  by_vars_text <- paste0(
    by_vars_text,
    # "The variables must be a unique key of the selected observations. \n \n",
    #

    "The expressions should be quoted using the `exprs()` function from ",
    "the rlang package.",
    "e.g. exprs(USUBJID, VISIT) \n \n"

  )

  if (!is.null(additional_dataset)) {
    by_vars_text <- paste0(
      by_vars_text,
      "The by variables must be a unique key of the selected observations in ",
      enumerate(additional_dataset), ". \n \n",
      "Variables from ", enumerate(additional_dataset),
      " can be renamed by naming the element, I.e. \n",
      "by_vars = exprs(<name in input dataset› = ‹name in additional dataset>), similar to the dplyr joins.\n \n"
    )
  }

  by_vars_text <- paste0(
    by_vars_text,
    "*Permitted Values*: list of variables created by exprs()"
    )

  return(by_vars_text)
}
