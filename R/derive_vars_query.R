#' Derive Query Variables
#'
#' @details This function can be used to derive CDISC variables such as
#'   `SMQzzNAM`, `SMQzzCD`, `SMQzzSC`, `SMQzzSCN`, and `CQzzNAM` in ADAE and
#'   ADMH, and variables such as `SDGzzNAM`, `SDGzzCD`, and `SDGzzSC` in ADCM.
#'   An example usage of this function can be found in the
#'   [OCCDS vignette](../articles/occds.html).
#'
#'   A query dataset is expected as an input to this function. See the
#'   [Queries Dataset Documentation vignette](../articles/queries_dataset.html)
#'   for descriptions, or call `data("queries")` for an example of a query dataset.
#'
#'   For each unique element in `VAR_PREFIX`, the corresponding "NAM"
#'   variable will be created. For each unique `VAR_PREFIX`, if `QUERY_ID` is
#'   not "" or NA, then the corresponding "CD" variable is created; similarly,
#'   if `QUERY_SCOPE` is not "" or NA, then the corresponding "SC" variable will
#'   be created; if `QUERY_SCOPE_NUM` is not "" or NA, then the corresponding
#'   "SCN" variable will be created.
#'
#'   For each record in `dataset`, the "NAM" variable takes the value of
#'   `QUERY_NAME` if the value of `TERM_NAME` or `TERM_ID` in `dataset_queries` matches
#'   the value of the respective TERM_LEVEL in `dataset`.
#'   Note that `TERM_NAME` in `dataset_queries` dataset may be NA only when `TERM_ID`
#'   is non-NA and vice versa.
#'   The "CD", "SC", and "SCN" variables are derived accordingly based on
#'   `QUERY_ID`, `QUERY_SCOPE`, and `QUERY_SCOPE_NUM` respectively,
#'   whenever not missing.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_queries A dataset containing required columns `VAR_PREFIX`,
#' `QUERY_NAME`, `TERM_LEVEL`, `TERM_NAME`, `TERM_ID`, and optional columns
#' `QUERY_ID`, `QUERY_SCOPE`, `QUERY_SCOPE_NUM`.
#'
#'   The content of the dataset will be verified by [assert_valid_queries()].
#'
#'   `create_query_data()` can be used to create the dataset.
#'
#' @author Ondrej Slama, Shimeng Huang
#'
#' @return The input dataset with query variables derived.
#'
#' @family der_occds
#' @keywords der_occds
#'
#' @seealso [create_query_data()] [assert_valid_queries()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' data("queries")
#' adae <- tribble(
#'   ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~AELLTCD,
#'   "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
#'   3, "Alanine aminotransferase abnormal", NA_character_, NA_integer_,
#'   "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
#'   5, "Basedow's disease", NA_character_, 1L,
#'   "03", "2020-06-07 23:59:59", "SOME TERM",
#'   2, "Some query", "Some term", NA_integer_,
#'   "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
#'   7, "Alveolar proteinosis", NA_character_, NA_integer_
#' )
#' derive_vars_query(adae, queries)
derive_vars_query <- function(dataset, dataset_queries) {
  assert_data_frame(dataset_queries)
  assert_valid_queries(dataset_queries, queries_name = deparse(substitute(dataset_queries)))
  assert_data_frame(dataset,
    required_vars = vars(!!!syms(unique(dataset_queries$TERM_LEVEL))),
    optional = FALSE
  )

  # replace all "" by NA
  dataset_queries <- dataset_queries %>%
    dplyr::mutate_if(is.character, function(x) {
      ifelse(x == "", NA_character_, x)
    })

  # names of new columns
  if ("QUERY_ID" %notin% names(dataset_queries)) {
    dataset_queries$QUERY_ID <- NA_integer_ # nolint
  }
  if ("QUERY_SCOPE" %notin% names(dataset_queries)) {
    dataset_queries$QUERY_SCOPE <- NA_integer_ # nolint
  }
  if ("QUERY_SCOPE_NUM" %notin% names(dataset_queries)) {
    dataset_queries$QUERY_SCOPE_NUM <- NA_integer_ # nolint
  }
  new_col_names <- dataset_queries %>%
    group_by(VAR_PREFIX) %>%
    mutate(
      NAM = paste0(VAR_PREFIX, "NAM"),
      CD = ifelse(!all(is.na(QUERY_ID)),
        paste0(VAR_PREFIX, "CD"), NA_character_
      ),
      SC = ifelse(!all(is.na(QUERY_SCOPE)),
        paste0(VAR_PREFIX, "SC"), NA_character_
      ),
      SCN = ifelse(!all(is.na(QUERY_SCOPE_NUM)),
        paste0(VAR_PREFIX, "SCN"), NA_character_
      )
    ) %>%
    ungroup() %>%
    select(NAM, CD, SC, SCN) %>%
    distinct() %>%
    pivot_longer(c(NAM, CD, SC, SCN), names_to = "key", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      order1 = stringr::str_extract(value, "^[a-zA-Z]{2,3}"),
      order2 = stringr::str_extract(value, "\\d{2}"),
      order3 = as.integer(factor(key, levels = c("NAM", "CD", "SC", "SCN")))
    ) %>%
    arrange(desc(order1), order2, order3) %>%
    pull(value)

  # queries restructured
  queries_wide <- dataset_queries %>%
    mutate(
      TERM_NAME = toupper(.data$TERM_NAME),
      VAR_PREFIX_NAM = paste0(.data$VAR_PREFIX, "NAM")
    ) %>%
    pivot_wider(names_from = .data$VAR_PREFIX_NAM, values_from = .data$QUERY_NAME) %>%
    mutate(VAR_PREFIX_CD = paste0(.data$VAR_PREFIX, "CD")) %>%
    pivot_wider(names_from = .data$VAR_PREFIX_CD, values_from = .data$QUERY_ID) %>%
    mutate(VAR_PREFIX_SC = paste0(.data$VAR_PREFIX, "SC")) %>%
    pivot_wider(names_from = .data$VAR_PREFIX_SC, values_from = .data$QUERY_SCOPE) %>%
    mutate(VAR_PREFIX_SCN = paste0(.data$VAR_PREFIX, "SCN")) %>%
    pivot_wider(names_from = .data$VAR_PREFIX_SCN, values_from = .data$QUERY_SCOPE_NUM) %>%
    select(-VAR_PREFIX) %>%
    # determine join column based on type of TERM_LEVEL
    # numeric -> TERM_ID, character -> TERM_NAME, otherwise -> error
    mutate(
      tmp_col_type = vapply(dataset[.data$TERM_LEVEL], typeof, character(1)),
      TERM_NAME_ID = case_when(
        .data$tmp_col_type == "character" ~ .data$TERM_NAME,
        .data$tmp_col_type %in% c("double", "integer") ~ as.character(.data$TERM_ID),
        TRUE ~ NA_character_
      )
    )

  # throw error if any type of column is not character or numeric
  if (any(is.na(queries_wide$TERM_NAME_ID))) {
    idx <- is.na(queries_wide$TERM_NAME_ID)
    dat_incorrect_type <- dataset[queries_wide$TERM_LEVEL[idx]]
    msg <- paste0(
      paste0(
        colnames(dat_incorrect_type),
        " is of type ",
        vapply(dat_incorrect_type, typeof, character(1)),
        collapse = ", "
      ),
      ", numeric or character is required"
    )
    abort(msg)
  }

  # prepare input dataset for joining
  static_cols <- setdiff(names(dataset), unique(dataset_queries$TERM_LEVEL))
  # if dataset does not have a unique key, create a temp one
  no_key <- dataset %>%
    select(!!!syms(static_cols)) %>%
    distinct()
  if (nrow(no_key) != nrow(dataset)) {
    dataset$temp_key <- seq_len(nrow(dataset))
    static_cols <- c(static_cols, "temp_key")
  }

  # Keep static variables - will add back on once non-static vars fixed
  df_static <- dataset %>% select(static_cols)

  # Change non-static numeric vars to character
  df_fix_numeric <- dataset %>%
    select(-static_cols) %>%
    mutate_if(is.numeric, as.character)


  joined <- cbind(df_static, df_fix_numeric) %>%
    pivot_longer(-static_cols, names_to = "TERM_LEVEL", values_to = "TERM_NAME_ID") %>%
    drop_na(.data$TERM_NAME_ID) %>%
    mutate(TERM_NAME_ID = toupper(.data$TERM_NAME_ID))

  # join restructured queries to input dataset
  joined <- joined %>%
    inner_join(queries_wide, by = c("TERM_LEVEL", "TERM_NAME_ID")) %>%
    select(!!!syms(c(static_cols, new_col_names))) %>%
    dplyr::group_by_at(static_cols) %>%
    dplyr::summarise_all(~ dplyr::first(na.omit(.))) %>%
    ungroup()

  # join queries to input dataset
  derive_vars_merged(dataset, dataset_add = joined, by_vars = vars(!!!syms(static_cols))) %>%
    select(-starts_with("temp_"))
}

#' Verify if a Dataset Has the Required Format as Queries Dataset.
#'
#' @details Check if the dataset has the following columns
#' - `VAR_PREFIX`, e.g., SMQ01, CQ12
#' - `QUERY_NAME`, non NA, must be unique per each `VAR_PREFIX`
#' - `QUERY_ID`, could be NA, must be unique per each `VAR_PREFIX`
#' - `QUERY_SCOPE`, 'BROAD', 'NARROW', or NA
#' - `QUERY_SCOPE_NUM`, 1, 2, or NA
#' - `TERM_LEVEL`, e.g., `"AEDECOD"`, `"AELLT"`, `"AELLTCD"`, ...
#' - `TERM_NAME`, character, could be NA only at those observations
#' where `TERM_ID` is non-NA
#' - `TERM_ID`, integer, could be NA only at those observations
#' where `TERM_NAME` is non-NA
#'
#' @param queries A data.frame.
#'
#' @param queries_name Name of the queries dataset, a string.
#'
#' @author Shimeng Huang, Ondrej Slama
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @export
#'
#' @return The function throws an error if any of the requirements not met.
#'
#' @examples
#' data("queries")
#' assert_valid_queries(queries, "queries")
assert_valid_queries <- function(queries, queries_name) {
  # check required columns
  assert_has_variables(
    queries,
    c("VAR_PREFIX", "QUERY_NAME", "TERM_LEVEL", "TERM_NAME", "TERM_ID")
  )

  # check duplicate rows
  signal_duplicate_records(queries, by_vars = quos(!!!syms(colnames(queries))))

  # check illegal prefix category
  is_good_prefix <- grepl("^[a-zA-Z]{2,3}", queries$VAR_PREFIX)
  if (!all(is_good_prefix)) {
    abort(
      paste0(
        "`VAR_PREFIX` in `", queries_name,
        "` must start with 2-3 letters.. Problem with ",
        enumerate(unique(queries$VAR_PREFIX[!is_good_prefix])),
        "."
      )
    )
  }

  # check illegal prefix number
  query_num <- sub("[[:alpha:]]+", "", queries$VAR_PREFIX)
  is_bad_num <- nchar(query_num) != 2 | is.na(as.numeric(query_num))
  if (any(is_bad_num)) {
    abort(
      paste0(
        "`VAR_PREFIX` in `", queries_name,
        "` must end with 2-digit numbers. Issue with ",
        enumerate(unique(queries$VAR_PREFIX[is_bad_num])),
        "."
      )
    )
  }

  # check illegal query name
  if (any(queries$QUERY_NAME == "") | any(is.na(queries$QUERY_NAME))) {
    abort(paste0(
      "`QUERY_NAME` in `", queries_name,
      "` cannot be empty string or NA."
    ))
  }

  # check query id is numeric
  if ("QUERY_ID" %in% names(queries) && !is.numeric(queries$QUERY_ID)) {
    abort(paste0(
      "`QUERY_ID` in `", queries_name,
      "` should be numeric."
    ))
  }

  # check illegal query scope
  if ("QUERY_SCOPE" %in% names(queries) &&
    any(unique(queries$QUERY_SCOPE) %notin% c("BROAD", "NARROW", "", NA_character_))) {
    abort(paste0(
      "`QUERY_SCOPE` in `", queries_name,
      "` can only be 'BROAD', 'NARROW' or `NA`."
    ))
  }

  # check illegal query scope number
  if ("QUERY_SCOPE_NUM" %in% names(queries)) {
    is_bad_scope_num <- queries$QUERY_SCOPE_NUM %notin% c(1, 2, NA_integer_)
    if (any(is_bad_scope_num)) {
      abort(
        paste0(
          "`QUERY_SCOPE_NUM` in `", queries_name,
          "` must be one of 1, 2, or NA. Issue with ",
          enumerate(unique(queries$QUERY_SCOPE_NUM[is_bad_scope_num])),
          "."
        )
      )
    }
  }

  # check illegal term name
  if (any(is.na(queries$TERM_NAME) & is.na(queries$TERM_ID)) |
    any(queries$TERM_NAME == "" & is.na(queries$TERM_ID))) {
    abort(paste0(
      "Either `TERM_NAME` or `TERM_ID` need to be specified",
      " in `", queries_name, "`. ",
      "They both cannot be NA or empty."
    ))
  }

  # each VAR_PREFIX must have unique QUERY_NAME, QUERY_ID if the columns exist
  count_unique <- queries %>%
    group_by(VAR_PREFIX) %>%
    dplyr::summarise(
      n_qnam = length(unique(QUERY_NAME)),
      n_qid = ifelse("QUERY_ID" %in% names(queries),
        length(unique(QUERY_ID)), 0
      )
    ) %>%
    ungroup()

  if (any(count_unique$n_qnam > 1)) {
    idx <- which(count_unique$n_qnam > 1)
    abort(paste0(
      "In `", queries_name, "`, `QUERY_NAME` of '",
      paste(count_unique$VAR_PREFIX[idx], collapse = ", "),
      "' is not unique."
    ))
  }

  if (any(count_unique$n_qid > 1)) {
    idx <- which(count_unique$n_qid > 1)
    abort(paste0(
      "In `", queries_name, "`, `QUERY_ID` of '",
      paste(count_unique$VAR_PREFIX[idx], collapse = ", "),
      "' is not unique."
    ))
  }
}
