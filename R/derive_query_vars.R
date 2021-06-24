#' Derive query variables.
#'
#' @details For each unique element in `VAR_PREFIX`, the corresponding "NAM"
#'   variable will be created. For each unique `VAR_PREFIX`, if `QUERY_ID` is
#'   not "" or NA, then the corresponding "CD" variable is created; similarly,
#'   if `QUERY_SCOPE` is not "" or NA, then the corresponding "SC" variable will
#'   be created; if `QUERY_SCOPE_NUM` is not "" or NA, then the corresponding
#'   "SCN" variable will be created.
#'
#'   For each record in `dataset`, the "NAM" variable takes the value of
#'   `QUERY_NAME` if the value of `TERM_NAME` in `queries` matches the value
#'   in the `TERM_LEVEL` column of `dataset`. The "CD", "SC", and "SCN"
#'   variables are derived accordingly based on `QUERY_ID`, `QUERY_SCOPE`, and
#'   `QUERY_SCOPE_NUM` respectively, whenever not missing.
#'
#' @param dataset Input dataset.
#'
#' @param queries A data.frame containing required columns `VAR_PREFIX`,
#' `QUERY_NAME`, `TERM_LEVEL`, `TERM_NAME`, and optional columns
#' `QUERY_ID`, `QUERY_SCOPE`, `QUERY_SCOPE_NUM`.
#'
#'   The content of the dataset will be verified by [assert_valid_queries()].
#'
#' @author Ondrej Slama, Shimeng Huang
#'
#' @return The input dataset with query variables derived.
#'
#' @keywords adae adcm derivations
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
#' 7, "Alveolar proteinosis", NA_character_)
#' derive_query_vars(adae, queries)
derive_query_vars <- function(dataset, queries) {

  assert_that(
    is.data.frame(dataset),
    is.data.frame(queries)
  )
  assert_valid_queries(queries, deparse(substitute(queries)))

  # replace all "" by NA
  queries <- queries %>%
    dplyr::mutate_if(is.character, function(x) {
      ifelse(x == "", NA_character_, x)})

  # names of new columns
  if ("QUERY_ID" %!in% names(queries)) {
    queries$QUERY_ID <- NA_integer_ # nolint
  }
  if ("QUERY_SCOPE" %!in% names(queries)) {
    queries$QUERY_SCOPE <- NA_integer_ # nolint
  }
  if ("QUERY_SCOPE_NUM" %!in% names(queries)) {
    queries$QUERY_SCOPE_NUM <- NA_integer_ # nolint
  }
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
  new_col_names <- c(paste0(unique(queries$VAR_PREFIX), "NAM"), new_col_names)

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
  # if dataset does not have a unique key, create a temp one
  no_key <- dataset %>% select(static_cols) %>% distinct() %>% nrow() != nrow(dataset)
  if (no_key) {
    dataset$temp_key <- seq_len(nrow(dataset))
    static_cols <- c(static_cols, "temp_key")
  }
  joined <- dataset %>%
    gather(key = "TERM_LEVEL", value = "TERM_NAME", -static_cols) %>%
    drop_na(.data$TERM_NAME)

  # join restructured queries to input dataset
  joined <- joined %>%
    mutate(TERM_NAME_UPPER = toupper(.data$TERM_NAME)) %>%
    dplyr::inner_join(queries_wide, by = c("TERM_LEVEL", "TERM_NAME_UPPER" = "TERM_NAME")) %>%
    select(static_cols, new_col_names) %>%
    dplyr::group_by_at(static_cols) %>%
    dplyr::summarise_all(~dplyr::first(na.omit(.)))

  # join queries to input dataset
  left_join(dataset, joined,
            by = static_cols) %>%
    select(-starts_with("temp_"))
}

#' Verify if a dataset has the required format as queries dataset.
#'
#' @details Check if the dataset has the following columns
#' - `VAR_PREFIX`, e.g., SMQ01, CQ12
#' - `QUERY_NAME`, non NULL, must be unique per each `VAR_PREFIX`
#' - `QUERY_ID`, could be NULL, must be unique per each `VAR_PREFIX`
#' - `QUERY_SCOPE`, 'BROAD', 'NARROW', or NULL
#' - `QUERY_SCOPE_NUM`, 1, 2, or NA
#' - `TERM_LEVEL`, e.g., AEDECOD, AELLT, ...
#' - `TERM_NAME`, non NULL
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
assert_valid_queries <- function(queries, queries_name) {

  # check required columns
  required_cols <-  c("VAR_PREFIX", "QUERY_NAME", "TERM_LEVEL", "TERM_NAME")
  is_missing <- required_cols %!in% names(queries)
  if (any(is_missing)) {
    missing_vars <- required_cols[is_missing]
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
  bad_prefix <- nchar(sub("[^[:alpha:]]+", "", queries$VAR_PREFIX)) > 3
  if (sum(bad_prefix) == 1L) {
    err_msg <- paste0("`VAR_PREFIX` in `", queries_name,
                      "` must start with 2-3 letters.. Problem with `", bad_prefix, "`.")
    abort(err_msg)
  }
  else if (sum(bad_prefix) > 1L) {
    err_msg <- paste0(
      "`VAR_PREFIX` in `", queries_name,
      "` must start with 2-3 letters.. Problem with ",
      enumerate(bad_prefix),
      ".")
    abort(err_msg)
  }

  # check illegal prefix number
  query_num <- sub("[[:alpha:]]+", "", queries$VAR_PREFIX)
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
  if ("QUERY_ID" %in% names(queries) && !is.numeric(queries$QUERY_ID)) {
    abort(paste0("`QUERY_ID` in `", queries_name,
                 "` should be numeric."))
  }

  # check illegal query scope
  if ("QUERY_SCOPE" %in% names(queries) &&
      any(unique(queries$QUERY_SCOPE) %!in% c("BROAD", "NARROW", "", NA_character_))) {
    abort(paste0("`QUERY_SCOPE` in `", queries_name,
                 "` can only be 'BROAD', 'NARROW' or `NA`."))
  }

  # check illegal query scope number
  if ("QUERY_SCOPE_NUM" %in% names(queries)) {
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
  }

  # check illegal term name
  if (any(queries$TERM_NAME == "") | any(is.na(queries$TERM_NAME))) {
    abort(paste0("`TERM_NAME` in `", queries_name,
                 "` cannot be empty string or NA."))
  }

  # each VAR_PREFIX must have unique QUERY_NAME, QUERY_ID if the columns exist
  count_unique <- queries %>%
    group_by(VAR_PREFIX) %>%
    dplyr::summarise(n_qnam = length(unique(QUERY_NAME)),
                     n_qid = ifelse("QUERY_ID" %in% names(queries),
                                    length(unique(QUERY_ID)), 0)) %>%
    ungroup()
  for (ii in seq_len(nrow(count_unique))) {
    if (count_unique[ii, ]$n_qnam > 1) {
      abort(paste0("In `", queries_name, "`, `QUERY_NAME` of '",
                   count_unique$VAR_PREFIX[ii], "' is not unique."))
    }
    if (count_unique[ii, ]$n_qid > 1) {
      abort(paste0("In `", queries_name, "`, `QUERY_ID` of '",
                   count_unique$VAR_PREFIX[ii], "' is not unique."))
    }
  }
}
