#' Derive query variables.
#'
#' @details For each unique element in `VAR_PREFIX`, the corresponding "NAM"
#'   variable is created. For each unique `VAR_PREFIX`, if `QUERY_ID` is not ""
#'   or NA, then the corresponding "CD" variable is created ; similarly, if
#'   `QUERY_SCOPE`is not "" or NA, then the corresponding "SC" variable is
#'   created.
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
#' @author Ondrej Slama, Shimeng Huang
#'
#' @keywords adae adcm
#'
#' @seealso [assert_valid_queries()]
#'
#' @export
#'
#' @examples
#' data("queries")
#' adae <- tibble::tribble(
#' ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT,
#' "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
#' 3, "Alanine aminotransferase abnormal", NA_character_,
#' "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
#' 5, "Basedow's disease", NA_character_,
#' "03", "2020-06-07 23:59:59", "SOME TERM",
#' 2, "Some query", "Some term",
#' "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
#' 7, "Alveolar proteinosis", NA_character_
#' )
#' derive_query_vars(adae, queries)
derive_query_vars <- function(dataset, queries) {

  assert_that(
    is.data.frame(dataset),
    is.data.frame(queries)
  )
  assert_valid_queries(queries, deparse(substitute(queries)))

  # replace all "" by NA
  queries <- queries %>%
    dplyr::mutate_if(is.character, function(x) {ifelse(x == "", NA_character_, x)})
    # mutate(across(where(is.character), ~na_if(., ""))) %>% # (not available for older dplyr)

  # names of new columns
  new_col_names <- queries %>%
    group_by(VAR_PREFIX) %>%
    mutate(CD = ifelse(!all(is.na(QUERY_ID)),
                       paste0(VAR_PREFIX, "CD"), NA_character_),
           SC = ifelse(!all(is.na(QUERY_SCOPE)),
                       paste0(VAR_PREFIX, "SC"), NA_character_),
           SCN = ifelse(!all(is.na(QUERY_SCOPE_NUM)),
                        paste0(VAR_PREFIX, "SCN"), NA_character_)) %>%
    ungroup() %>%
    select(CD, SC, SCN) %>%
    gather() %>%
    filter(!is.na(value)) %>%
    pull(value) %>%
    unique()
  new_cols_names <- c(paste0(unique(queries$VAR_PREFIX), "NAM"), new_col_names)

  # queries restructured
  queries_wide <- queries %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME),
           VAR_PREFIX_NAM = paste0(.data$VAR_PREFIX, "NAM")) %>%
    spread(.data$VAR_PREFIX_NAM, .data$QUERY_NAME) %>%
    mutate(VAR_PREFIX_CD = paste0(.data$VAR_PREFIX, "CD")) %>%
    spread(.data$VAR_PREFIX_CD, .data$QUERY_ID) %>%
    mutate(VAR_PREFIX_SC = paste0(.data$VAR_PREFIX, "SC")) %>%
    spread(.data$VAR_PREFIX_SC, .data$QUERY_SCOPE)  %>%
    mutate(VAR_PREFIX_SCN = paste0(.data$VAR_PREFIX, "SCN")) %>%
    spread(.data$VAR_PREFIX_SCN, .data$QUERY_SCOPE_NUM)

  queries_wide <- queries_wide %>%
    select(-VAR_PREFIX) %>%
    mutate(TERM_NAME = toupper(.data$TERM_NAME))

  # prepare input dataset for joining
  static_cols <- setdiff(names(dataset), unique(queries$TERM_LEVEL))
  joined <- dataset %>%
    gather(key = "TERM_LEVEL", value = "TERM_NAME", -static_cols) %>%
    drop_na(.data$TERM_NAME)
  term_cols <- unique(joined$TERM_LEVEL)

  # join restructured queries to input dataset
  joined <- joined %>%
    mutate(TERM_NAME_UPPER = toupper(.data$TERM_NAME)) %>%
    dplyr::inner_join(queries_wide, by = c("TERM_LEVEL", "TERM_NAME_UPPER" = "TERM_NAME"))

  # join queries to input dataset
  left_join(dataset, select(joined, static_cols, new_cols_names), by = static_cols)
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

  # check required columns
  is_missing <- c("VAR_PREFIX", "QUERY_NAME",
                  "QUERY_ID", "QUERY_SCOPE", "QUERY_SCOPE_NUM",
                  "TERM_LEVEL", "TERM_NAME") %!in% names(queries)
  if (any(is_missing)) {
    missing_vars <- names(queries)[is_missing]
    if (length(missing_vars) == 1L) {
      err_msg <- paste0("Required variable in `", missing_vars,
                        "` is missing in `", queries_name, "`.")
    } else {
      err_msg <- paste0(
        "Required variables ",
        enumerate(missing_vars),
        " are missing in `", queries_name, "`."
      )
    }
    abort(err_msg)
  }

  # check duplicate rows
  if (nrow(queries) != nrow(queries %>% dplyr::distinct())) {
    abort("`", queries_name, "` should not have duplicate rows.")
  }

  # check illegal prefix category
  bad_prefix <- !grepl("^(SMQ|CQ|SDQ)", queries$VAR_PREFIX)
  if (any(bad_prefix)) {
    abort(paste0("`VAR_PREFIX` in `", queries_name,
                 "` must start with one of 'SMQ', 'SDQ', or 'CQ'."))
    if (length(bad_prefix) == 1L) {
      err_msg <- paste0("`VAR_PREFIX` in `", queries_name,
                        "` must start with one of 'SMQ', 'SDQ', or 'CQ'. Problem with `", bad_prefix, "`.")
    } else {
      err_msg <- paste0(
        "`VAR_PREFIX` in `", queries_name,
        "` must start with one of 'SMQ', 'SDQ', or 'CQ'. Problem with ",
        enumerate(bad_prefix),
        ".")
    }
    abort(err_msg)
  }

  # check illegal prefix number
  query_num <- gsub("^(SMQ|CQ|SDQ)", "", queries$VAR_PREFIX)
  is_bad_num <- nchar(query_num) != 2 | is.na(as.numeric(query_num))
  if (any(is_bad_num)) {
    bad_nums <- unique(queries$VAR_PREFIX[is_bad_num])
    if (length(bad_nums) == 1L) {
      err_msg <- paste0("`VAR_PREFIX` in `", queries_name,
                        "` must end with 2-digit numbers. Issue with `", bad_nums, "`.")
    } else {
      err_msg <- paste0(
        "`VAR_PREFIX` in `", queries_name,
        "` must end with 2-digit numbers. Issue with ",
        enumerate(bad_nums),
        ".")
    }
    abort(err_msg)
  }

  # check illegal query name
  if (any(queries$QUERY_NAME == "") | any(is.na(queries$QUERY_NAME))) {
    abort(paste0("`QUERY_NAME` in `", queries_name,
                 "` cannot be empty string or NA."))
  }

  # check query id is numeric
  if (class(queries$QUERY_ID) != "numeric") {
    abort(paste0("`QUERY_ID` in `", queries_name,
                 "` should be numeric."))
  }

  # check illegal query scope
  if (any(unique(queries$QUERY_SCOPE) %!in% c("BROAD", "NARROW", "", NA_character_))) {
    abort(paste0("`QUERY_SCOPE` in `", queries_name,
                 "` can only be 'BROAD', 'NARROW' or `NA`."))
  }

  # check illegal query scope number
  is_bad_scope_num <- queries$QUERY_SCOPE_NUM %!in% c(1, 2, NA_integer_)
  if (any(is_bad_scope_num)) {
    bad_scope_nums <- unique(queries$QUERY_SCOPE_NUM[is_bad_scope_num])
    if (length(bad_scope_nums) == 1L) {
      err_msg <- paste0("`QUERY_SCOPE_NUM` in `", queries_name,
                        "` must be one of 1, 2, or NA. Issue with `", bad_scope_nums, "`.")
    } else {
      err_msg <- paste0(
        "`QUERY_SCOPE_NUM` in `", queries_name,
        "` must be one of 1, 2, or NA. Issue with ",
        enumerate(bad_scope_nums),
        ".")
    }
    abort(err_msg)
  }

  # check illegal term name
  if (any(queries$TERM_NAME == "") | any(is.na(queries$TERM_NAME))) {
    abort(paste0("`TERM_NAME` in `", queries_name,
                 "` cannot be empty string or NA."))
  }

  # each VAR_PREFIX must have unique QUERY_NAME, QUERY_ID, and QUERY_SCOPE
  count_unique <- queries %>%
    group_by(VAR_PREFIX) %>%
    dplyr::summarise(n_qnam = length(unique(QUERY_NAME)),
                     n_qid = length(unique(QUERY_ID)),
                     n_qsc = length(unique(QUERY_SCOPE)))
  for (ii in 1:nrow(count_unique)) {
    if (count_unique[ii,]$n_qnam > 1) {
      abort(paste0("In `", queries_name, "`, `QUERY_NAME` of '",
                   count_unique$VAR_PREFIX[ii], "' is not unique."))
    }
    if (count_unique[ii,]$n_qid > 1) {
      abort(paste0("In `", queries_name, "`, `QUERY_ID` of '",
                   count_unique$VAR_PREFIX[ii], "' is not unique."))
    }
    if (count_unique[ii,]$n_qsc > 1) {
      abort(paste0("In `", queries_name, "`, `QUERY_SCOPE` of '",
                   count_unique$VAR_PREFIX[ii], "' is not unique."))
    }
  }
}

