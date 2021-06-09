#' Derive query variables.
#'
#' @details If `VAR_PREFIX` starts with "SMQ" or "SDG", then the corresponding
#'   "NAM", "CD", and "SC" variables are derived; if starts with "CQ", only
#'   "NAM" variable is derived.
#'
#' @param dataset Input dataset.
#'
#'   The columns specified by `dataset_keys` are expected.
#'
#' @param queries A data.frame containing the following required columns:
#' - VAR_PREFIX, e.g., SMQ01, CQ12
#' - QUERY_NAME, non NULL
#' - QUERY_ID, could be NULL
#' - QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
#' - TERM_LEVEL, e.g., AEDECOD, AELLT, ...
#' - TERM_NAME, non NULL
#'   This will be check by [validate_queries()].
#'
#' @param dataset_keys A vector of column names that can be used to identify
#'   a unique record in the `dataset`. e.g. for ADAE, this would be
#'   `c("USUBJID", "ASTDTM", "AETERM", "AESEQ")`.
#'
#' @author Shimeng Huang
#'
#' @keywords adae adcm
#'
#' @export
derive_query_vars <- function(dataset, queries, dataset_keys) {

  assert_has_variables(dataset, dataset_keys)
  assert_valid_queries(queries)

  # names of new columns
  new_cols_names <- sapply(queries$VAR_PREFIX, function(x) {
    if (grepl("SMQ", x) | grepl("SGD", x)) return(paste0(x, c("NAM", "CD", "SC")))
    else if (grepl("CQ", x)) return(paste0(x, "NAM"))
    else return("")
  }, USE.NAMES = FALSE)
  new_cols_names <- unlist(new_cols_names)

  # add in empty new columns in dataset
  dataset[new_cols_names] <- NA_character_

  # split by each level of query since the joining variable would be different
  queries_list <- split(queries, queries$TERM_LEVEL)
  for (queries_one_level in queries_list) {
    term_level <- unique(queries_one_level$TERM_LEVEL) # the column name to be joined by in ADAE
    dataset_temp <- dataset %>%
      select(dataset_keys, term_level)

    # TODO: better way to match ignoring case and not use fuzzyjoin? see below in left_join too
    dataset_temp[[paste0(term_level, "_temp")]] <- toupper(dataset_temp[[term_level]])
    dataset_temp <- dataset_temp %>%
      left_join(queries_one_level %>%
                  mutate(TERM_NAME = toupper(TERM_NAME)),
                by = setNames("TERM_NAME", paste0(term_level, "_temp"))) %>%
      filter(!is.na(VAR_PREFIX)) # remove unmatched ones for this term_level
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
      select(dataset_keys, new_cols_names, term_level)

    dataset <- dataset %>%
      anti_join(dataset_temp, by = dataset_keys) %>%
      bind_rows(dataset_temp)
  }
  dataset
}

#' Verify if a dataset has the required format as queries dataset.
#'
#' @details Check if the dataset has the following columns
#' - VAR_PREFIX, e.g., SMQ01, CQ12
#' - QUERY_NAME, non NULL
#' - QUERY_ID, could be NULL
#' - QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
#' - TERM_LEVEL, e.g., AEDECOD, AELLT, ...
#' - TERM_NAME, non NULL
#'   This will be check by [validate_queries()].
#'
#' @param A data.frame.
#'
#' @author Shimeng Huang
#'
#' @export
#'
#' @return The function throws an error if any of the requirements not met.
assert_valid_queries <- function(queries) {
  if (any(c("VAR_PREFIX", "QUERY_NAME", "QUERY_ID", "QUERY_SCOPE",
            "TERM_LEVEL", "TERM_NAME") %!in% names(queries))) {
    abort("Missing required column(s) in `queries`.")
  }

  if (length(unique(queries$VAR_PREFIX)) != nrow(queries)) {
    abort("`VAR_PREFIX` in `queries` cannot have duplicates.")
  }

  # TODO: not sure if this one needed?
  query_num <- substr(queries$VAR_PREFIX,
                      start = nchar(queries$VAR_PREFIX)-1,
                      stop = nchar(queries$VAR_PREFIX))
  if (anyNA(as.numeric(query_num))) {
    abort("`VAR_PREFIX` in `queries` must end with 2-digit numbers.")
  }

  query_prefix <- unique(gsub('.{2}$', '', queries$VAR_PREFIX))
  if (any(query_prefix %!in% c("SMQ", "SDQ", "CQ"))) {
    abort("`VAR_PREFIX` in `queries` must start with one of 'SMQ', 'SDQ', or 'CQ'.")
  }

  if (any(queries$QUERY_NAME == "") | any(is.na(queries$QUERY_NAME))) {
    abort("`QUERY_NAME` in `queries` cannot be empty string or NA.")
  }

  # TODO: allow only NA instead of empty strings?
  if (any(unique(queries$QUERY_SCOPE) %!in% c("BROAD", "NARROW", NA_character_))) {
    abort("`QUERY_SCOPE` in `queries` can only be 'BROAD', 'NARROW' or `NA_character_`.")
  }

  if (any(queries$TERM_NAME == "") | any(is.na(queries$TERM_NAME))) {
    abort("`TERM_NAME` in `queries` cannot be empty string or NA.")
  }
}

