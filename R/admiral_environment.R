admiral_environment <- new.env(parent = emptyenv())

# admiral_options.R ----
## set_admiral_options
admiral_environment$admiral_options <- list(
  # future_input = vars(...), nolint
  subject_keys = vars(STUDYID, USUBJID)
)

# To enhance features and add inputs as necessary 

# 1. Add additional options such as future_input as shown commented above
# 2. Update @params with future_input in set_admiral_options roxygen documentation
# 3. Add future_input into set_admiral_options() formals and body

# derive_merged.R ----
## derive_vars_merged_lookup
admiral_environment$nmap <- NULL

# duplicates.R ----
## signal_duplicate_records
admiral_environment$duplicates <- NULL
