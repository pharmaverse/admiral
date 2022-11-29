#' Environment Objects
#'
#' @details
#' Once in a while, we may encounter "locked binding for 'xxx'." errors
#' during the development process while building out functions. This may arise because
#' we want to create dynamic data/objects based on user-inputs that need modification
#' at points in time after the package has been loaded. To manage such data or objects,
#' R has a data structure known as an 'environment'. These environment objects are created
#' at build time, but can be populated with values after the package has been loaded and
#' update those values over the course of an R session. For more details how environments work,
#' see relevant sections on environments in R Packages and Advanced R textbooks for more details.
#' @noRd

admiral_environment <- new.env(parent = emptyenv())

# See respective ...R page for usage

# admiral_options.R ----
## set_admiral_options
admiral_environment$admiral_options <- list(
  # future_input = vars(...), nolint
  subject_keys = vars(STUDYID, USUBJID)
)

# To enhance features and add inputs as necessary

# 1. Add additional options such as future_input as shown in comment above
# 2. Update @params with future_input in set_admiral_options roxygen documentation
# 3. Add future_input into set_admiral_options() formals and body

# derive_merged.R ----
## derive_vars_merged_lookup
admiral_environment$nmap <- NULL

# duplicates.R ----
## signal_duplicate_records
admiral_environment$duplicates <- NULL
