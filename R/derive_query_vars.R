#' Derive query variables.
#'
#' @details If `VAR_PREFIX` starts with "SMQ" or "SDG", then the corresponding
#'   "NAM", "CD", and "SC" variables are derived; if starts with "CQ", only
#'   "NAM" variable is derived.
#'
#' @param dataset Input dataset.
#' @param queries A data.frame containing the following required columns:
#' - VAR_PREFIX, e.g., SMQ01, CQ12
#' - QUERY_NAME, non NULL
#' - QUERY_ID, could be NULL
#' - QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
#' - TERM_LEVEL, e.g., AEDECOD, AELLT, ...
#' - TERM_NAME, non NULL
#'
#' @keywords adae adcm
#'
#' @export
derive_query_vars <- function(dataset, queries) {

  adae_keys <- c("USUBJID", "ASTDTM", "AETERM", "AESEQ")
  assert_has_variables(dataset, adae_keys)

  # TODO: also need to check no duplicates in VAR_PREFIX, no NULL in certain columns etc.
  assert_has_variables(queries,
                       c("VAR_PREFIX", "QUERY_NAME", "QUERY_ID", "QUERY_SCOPE",
                         "TERM_LEVEL", "TERM_NAME"))


  # names of new columns
  new_cols_names <- sapply(queries$VAR_PREFIX, function(x) {
    if (grepl("SMQ", x) | grepl("SGD", x)) return(paste0(x, c("NAM", "CD", "SC")))
    else if (grepl("CQ", x)) return(paste0(x, "NAM"))
    else return("")
  }, USE.NAMES = FALSE)
  new_cols_names <- unlist(new_cols_names)
  dataset[new_cols_names] <- NA_character_

  # split by each level of query since the joining variables are different
  queries_list <- split(queries, queries$TERM_LEVEL)
  for (queries_one_level in queries_list) {
    term_level <- unique(queries_one_level$TERM_LEVEL) # the column name to be joined by in ADAE
    dataset_temp <- dataset %>%
      select(adae_keys, term_level)

    # TODO: better way to match ignoring case and not use fuzzyjoin? see below in left_join too
    dataset_temp[[paste0(term_level, "_temp")]] <- toupper(dataset_temp[[term_level]])
    dataset_temp <- dataset_temp %>%
      left_join(queries_one_level %>%
                  mutate(TERM_NAME = toupper(TERM_NAME)),
                by = setNames("TERM_NAME", paste0(term_level, "_temp")))
    dataset_temp[new_cols_names] <- NA_character_

    # go through each row and assign values to the correct columns
    # TODO: better ways to do this?
    for (ii in 1:nrow(dataset_temp)) {
      var_prefix <- dataset_temp$VAR_PREFIX[ii]
      dataset_temp[[paste0(var_prefix, "NAM")]][ii] <- dataset_temp$QUERY_NAME[ii]
      if (grepl("SMQ", var_prefix) | grepl("SGD", var_prefix) ) {
        dataset_temp[[paste0(var_prefix, "SC")]][ii] <- dataset_temp$QUERY_SCOPE[ii]
        dataset_temp[[paste0(var_prefix, "CD")]][ii] <- dataset_temp$QUERY_ID[ii]
      }
    }

    dataset_temp <- dataset_temp %>%
      select(adae_keys, new_cols_names, term_level)

    dataset <- dataset %>%
      anti_join(dataset_temp, by = adae_keys) %>%
      bind_rows(dataset_temp)
  }
  dataset
}

