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
#'   This will be check by [assert_valid_queries()].
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
  new_cols_names <- lapply(queries$VAR_PREFIX, function(x) {
    if (grepl("SMQ", x) | grepl("SGD", x)) return(paste0(x, c("NAM", "CD", "SC")))
    else if (grepl("CQ", x)) return(paste0(x, "NAM"))
    else return("")
  })
  new_cols_names <- unlist(new_cols_names)
  # Note: in case something not matched? should this be checked in `queries`
  new_cols_names <- new_cols_names[new_cols_names != ""]

  # queries restructured
  qrs_tmp <- queries %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME),
           VAR_PREFIX_NAM = paste0(.data$VAR_PREFIX, "NAM")) %>%
    spread(.data$VAR_PREFIX_NAM, .data$QUERY_NAME) %>%
    mutate(VAR_PREFIX_CD = ifelse(grepl("^CQ", .data$VAR_PREFIX),
                                  "tmp_drop",
                                  paste0(.data$VAR_PREFIX, "CD"))) %>%
    spread(.data$VAR_PREFIX_CD, .data$QUERY_ID) %>%
    mutate(VAR_PREFIX_SC = ifelse(grepl("^CQ", .data$VAR_PREFIX),
                                  "tmp_drop2",
                                  paste0(.data$VAR_PREFIX, "SC"))) %>%
    spread(.data$VAR_PREFIX_SC, .data$QUERY_SCOPE) %>%
    select(-dplyr::matches("^tmp_drop"), -.data$VAR_PREFIX) %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME))

  # prepare input dataset for joining
  out <- dataset %>%
    gather(key = "TERM_LEVEL", value = "TERM_NAME", -dataset_keys) %>%
    drop_na(.data$TERM_NAME)
  term_cols <- unique(out$TERM_LEVEL)

  # join restructured queries to input dataset
  out <- out %>%
    mutate(TERM_NAME_UPPER = toupper(.data$TERM_NAME)) %>%
    dplyr::inner_join(qrs_tmp, by = c("TERM_LEVEL", "TERM_NAME_UPPER" = "TERM_NAME")) %>%
    spread(.data$TERM_LEVEL, .data$TERM_NAME) %>%
    select(dataset_keys, term_cols, new_cols_names, -.data$TERM_NAME_UPPER)

  out
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
#'
#' @param queries data.frame.
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

