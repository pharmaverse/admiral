#' Derive query variables.
#'
#' @details For each unique element in `VAR_PREFIX`, the corresponding "NAM"
#'   variable is created. For "SMQ" or "SGD", if `QUERY_ID` is not "" or NA,
#'   then the corresponding "CD" variable is created ; similarly, if `QUERY_SCOPE`
#'   is not "" or NA, then the corresponding "SC" variable is created.
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
#' @author Ondrej Slama, Shimeng Huang
#'
#' @keywords adae adcm
#'
#' @export
derive_query_vars <- function(dataset, queries, dataset_keys) {

  assert_has_variables(dataset, dataset_keys)
  assert_valid_queries(queries, deparse(substitute(queries)))

  # replace all "" by NA
  queries <- queries %>%
    dplyr::mutate_if(is.character, function(x) {ifelse(x == "", NA_character_, x)})
    # mutate(across(where(is.character), ~na_if(., ""))) %>% # (not available for older dplyr)

  # names of new columns
  # Note: currenlt as long as one of QUERY_ID in the group is not NA then
  #   "CD" variable is created, same for QUERY_SCOPE
  nam_names <- paste0(unique(queries$VAR_PREFIX), "NAM")
  cd_names <- queries %>%
    group_by(VAR_PREFIX) %>%
    filter(!all(is.na(QUERY_ID)) & !grepl("^CQ", VAR_PREFIX)) %>%
    ungroup() %>%
    pull(VAR_PREFIX) %>%
    unique()
  cd_names <- if_non_len0(cd_names, paste0(cd_names, "CD"))
  sc_names <- queries %>%
    group_by(VAR_PREFIX) %>%
    filter(!all(is.na(QUERY_SCOPE)) & !grepl("^CQ", VAR_PREFIX)) %>%
    ungroup() %>%
    pull(VAR_PREFIX) %>%
    unique()
  sc_names <- if_non_len0(sc_names, paste0(sc_names, "SC"))
  new_cols_names <- c(nam_names, cd_names, sc_names)

  # queries restructured
  queries_wide <- queries %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME),
           VAR_PREFIX_NAM = paste0(.data$VAR_PREFIX, "NAM")) %>%
    spread(.data$VAR_PREFIX_NAM, .data$QUERY_NAME)

  if (any(!is.na(queries$QUERY_ID))) {
    queries_wide <- queries_wide %>%
      mutate(VAR_PREFIX_CD = ifelse(grepl("^CQ", .data$VAR_PREFIX),
                                    "tmp_drop",
                                    paste0(.data$VAR_PREFIX, "CD"))) %>%
      spread(.data$VAR_PREFIX_CD, .data$QUERY_ID)
  }

  if (any(!is.na(queries$QUERY_SCOPE))) {
    queries_wide <- queries_wide %>%
      mutate(VAR_PREFIX_SC = ifelse(grepl("^CQ", .data$VAR_PREFIX),
                                    "tmp_drop2",
                                    paste0(.data$VAR_PREFIX, "SC"))) %>%
      spread(.data$VAR_PREFIX_SC, .data$QUERY_SCOPE)
  }

  queries_wide <- queries_wide %>%
    select(-dplyr::matches("^tmp_drop"), -.data$VAR_PREFIX) %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME))

  # prepare input dataset for joining
  joined <- dataset %>%
    gather(key = "TERM_LEVEL", value = "TERM_NAME", -dataset_keys) %>%
    drop_na(.data$TERM_NAME)
  term_cols <- unique(joined$TERM_LEVEL)

  # join restructured queries to input dataset
  joined <- joined %>%
    mutate(TERM_NAME_UPPER = toupper(.data$TERM_NAME)) %>%
    dplyr::inner_join(queries_wide, by = c("TERM_LEVEL", "TERM_NAME_UPPER" = "TERM_NAME"))

  # join queries to input dataset
  left_join(dataset, select(joined, dataset_keys, new_cols_names), by = dataset_keys)
}

#' Verify if a dataset has the required format as queries dataset.
#'
#' @details Check if the dataset has the following columns
#' - VAR_PREFIX, e.g., SMQ01, CQ12
#' - QUERY_NAME, non NULL, must be unique per each VAR_PREFIX
#' - QUERY_ID, could be NULL, must be unique per each VAR_PREFIX
#' - QUERY_SCOPE, 'BROAD', 'NARROW', or NULL, must be unique per each VAR_PREFIX
#' - TERM_LEVEL, e.g., AEDECOD, AELLT, ...
#' - TERM_NAME, non NULL
#'
#' @param queries A data.frame.
#'
#' @param queries_name Name of the queries dataset, a string.
#'
#' @author Shimeng Huang
#'
#' @export
#'
#' @return The function throws an error if any of the requirements not met.
assert_valid_queries <- function(queries,
                                 queries_name = deparse(substitute(queries))) {
  if (any(c("VAR_PREFIX", "QUERY_NAME",
            "QUERY_ID", "QUERY_SCOPE",
            "TERM_LEVEL", "TERM_NAME") %!in% names(queries))) {
    abort("Missing required column(s) in `queries`.")
  }

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

  if (any(unique(queries$QUERY_SCOPE) %!in% c("BROAD", "NARROW", "", NA_character_))) {
    abort("`QUERY_SCOPE` in `queries` can only be 'BROAD', 'NARROW' or `NA_character_`.")
  }

  if (any(queries$TERM_NAME == "") | any(is.na(queries$TERM_NAME))) {
    abort("`TERM_NAME` in `queries` cannot be empty string or NA.")
  }

  # each VAR_PREFIX must have unique QUERY_NAME, QUERY_ID, and QUERY_SCOPE
  count_unique <- queries %>%
    group_by(VAR_PREFIX) %>%
    dplyr::summarise(n_qnam = length(unique(QUERY_NAME)),
                     n_qid = length(unique(QUERY_ID)),
                     n_qsc = length(unique(QUERY_SCOPE)))
  for (ii in 1:nrow(count_unique)) {
    if (count_unique[ii,]$n_qnam > 1) {
      abort(paste0("QUERY_NAME of ", count_unique$VAR_PREFIX[ii], " is not unique."))
    }
    if (count_unique[ii,]$n_qid > 1) {
      abort(paste0("QUERY_ID of ", count_unique$VAR_PREFIX[ii], " is not unique."))
    }
    if (count_unique[ii,]$n_qsc > 1) {
      abort(paste0("QUERY_SCOPE of ", count_unique$VAR_PREFIX[ii], " is not unique."))
    }
  }
}

