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
