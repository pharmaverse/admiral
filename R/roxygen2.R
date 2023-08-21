# THIS FILE INCLUDES FUNCTION THAT ARE RE-USED THROUGHOUT
# roxygen2 COMMENTS TO GENERATE DOCUMENTATION TEXT

roxygen_param_dataset <- function(expected_vars_args) {
  paste0("Input dataset \n \n",
         "The variables specified by the ",
         enumerate(expected_vars_args),
         " argument(s) to be expected.")

}
