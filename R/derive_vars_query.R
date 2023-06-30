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
#'   For each unique element in `PREFIX`, the corresponding "NAM"
#'   variable will be created. For each unique `PREFIX`, if `GRPID` is
#'   not "" or NA, then the corresponding "CD" variable is created; similarly,
#'   if `SCOPE` is not "" or NA, then the corresponding "SC" variable will
#'   be created; if `SCOPEN` is not "" or NA, then the corresponding
#'   "SCN" variable will be created.
#'
#'   For each record in `dataset`, the "NAM" variable takes the value of
#'   `GRPNAME` if the value of `TERMNAME` or `TERMID` in `dataset_queries` matches
#'   the value of the respective SRCVAR in `dataset`.
#'   Note that `TERMNAME` in `dataset_queries` dataset may be NA only when `TERMID`
#'   is non-NA and vice versa.
#'   The "CD", "SC", and "SCN" variables are derived accordingly based on
#'   `GRPID`, `SCOPE`, and `SCOPEN` respectively,
#'   whenever not missing.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_queries A dataset containing required columns `PREFIX`,
#' `GRPNAME`, `SRCVAR`, `TERMNAME`, `TERMID`, and optional columns
#' `GRPID`, `SCOPE`, `SCOPEN`.
#'
#'   The content of the dataset will be verified by [assert_valid_queries()].
#'
#'   `create_query_data()` can be used to create the dataset.
#'
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
    required_vars = exprs(!!!syms(unique(dataset_queries$SRCVAR))),
    optional = FALSE
  )

  dataset_queries <- convert_blanks_to_na(dataset_queries)

  # names of new columns
  if ("GRPID" %notin% names(dataset_queries)) {
    dataset_queries$GRPID <- NA_integer_ # nolint
  }
  if ("SCOPE" %notin% names(dataset_queries)) {
    dataset_queries$SCOPE <- NA_integer_ # nolint
  }
  if ("SCOPEN" %notin% names(dataset_queries)) {
    dataset_queries$SCOPEN <- NA_integer_ # nolint
  }
  new_col_names <- dataset_queries %>%
    group_by(PREFIX) %>%
    mutate(
      NAM = paste0(PREFIX, "NAM"),
      CD = ifelse(!all(is.na(GRPID)),
        paste0(PREFIX, "CD"), NA_character_
      ),
      SC = ifelse(!all(is.na(SCOPE)),
        paste0(PREFIX, "SC"), NA_character_
      ),
      SCN = ifelse(!all(is.na(SCOPEN)),
        paste0(PREFIX, "SCN"), NA_character_
      )
    ) %>%
    ungroup() %>%
    select(NAM, CD, SC, SCN) %>%
    distinct() %>%
    pivot_longer(c(NAM, CD, SC, SCN), names_to = "key", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      order1 = str_extract(value, "^[a-zA-Z]{2,3}"),
      order2 = str_extract(value, "\\d{2}"),
      order3 = as.integer(factor(key, levels = c("NAM", "CD", "SC", "SCN")))
    ) %>%
    arrange(desc(order1), order2, order3) %>%
    pull(value)

  # queries restructured
  queries_wide <- dataset_queries %>%
    mutate(
      TERMNAME = toupper(TERMNAME),
      PREFIX_NAM = paste0(PREFIX, "NAM")
    ) %>%
    pivot_wider(names_from = PREFIX_NAM, values_from = GRPNAME) %>%
    mutate(PREFIX_CD = paste0(PREFIX, "CD")) %>%
    pivot_wider(names_from = PREFIX_CD, values_from = GRPID) %>%
    mutate(PREFIX_SC = paste0(PREFIX, "SC")) %>%
    pivot_wider(names_from = PREFIX_SC, values_from = SCOPE) %>%
    mutate(PREFIX_SCN = paste0(PREFIX, "SCN")) %>%
    pivot_wider(names_from = PREFIX_SCN, values_from = SCOPEN) %>%
    select(-PREFIX) %>%
    # determine join column based on type of SRCVAR
    # numeric -> TERMID, character -> TERMNAME, otherwise -> error
    mutate(
      tmp_col_type = vapply(dataset[SRCVAR], typeof, character(1)),
      TERM_NAME_ID = case_when(
        tmp_col_type == "character" ~ TERMNAME,
        tmp_col_type %in% c("double", "integer") ~ as.character(TERMID),
        TRUE ~ NA_character_
      )
    )

  # throw error if any type of column is not character or numeric
  if (any(is.na(queries_wide$TERM_NAME_ID))) {
    idx <- is.na(queries_wide$TERM_NAME_ID)
    dat_incorrect_type <- dataset[queries_wide$SRCVAR[idx]]
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
  static_cols <- setdiff(names(dataset), unique(dataset_queries$SRCVAR))
  # if dataset does not have a unique key, create a temp one
  no_key <- dataset %>%
    select(all_of(static_cols)) %>%
    distinct()
  if (nrow(no_key) != nrow(dataset)) {
    dataset$temp_key <- seq_len(nrow(dataset))
    static_cols <- c(static_cols, "temp_key")
  }

  # Keep static variables - will add back on once non-static vars fixed
  df_static <- dataset %>% select(all_of(static_cols))

  # Change non-static numeric vars to character
  df_fix_numeric <- dataset %>%
    select(-all_of(static_cols)) %>%
    mutate(across(where(is.numeric), as.character))


  joined <- cbind(df_static, df_fix_numeric) %>%
    pivot_longer(-all_of(static_cols), names_to = "SRCVAR", values_to = "TERM_NAME_ID") %>%
    drop_na(TERM_NAME_ID) %>%
    mutate(TERM_NAME_ID = toupper(TERM_NAME_ID))

  # join restructured queries to input dataset
  joined <- joined %>%
    inner_join(queries_wide, by = c("SRCVAR", "TERM_NAME_ID")) %>%
    select(!!!syms(c(static_cols, new_col_names))) %>%
    group_by_at(static_cols) %>%
    summarise_all(~ first(na.omit(.))) %>%
    ungroup()

  # join queries to input dataset
  derive_vars_merged(dataset, dataset_add = joined, by_vars = exprs(!!!syms(static_cols))) %>%
    select(-starts_with("temp_"))
}

#' Verify if a Dataset Has the Required Format as Queries Dataset.
#'
#' @details Check if the dataset has the following columns
#' - `PREFIX`, e.g., SMQ01, CQ12
#' - `GRPNAME`, non NA, must be unique per each `PREFIX`
#' - `GRPID`, could be NA, must be unique per each `PREFIX`
#' - `SCOPE`, 'BROAD', 'NARROW', or NA
#' - `SCOPEN`, 1, 2, or NA
#' - `SRCVAR`, e.g., `"AEDECOD"`, `"AELLT"`, `"AELLTCD"`, ...
#' - `TERMNAME`, character, could be NA only at those observations
#' where `TERMID` is non-NA
#' - `TERMID`, integer, could be NA only at those observations
#' where `TERMNAME` is non-NA
#'
#' @param queries A data.frame.
#'
#' @param queries_name Name of the queries dataset, a string.
#'
#'
#' @keywords other_advanced
#' @family other_advanced
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
  assert_data_frame(
    queries,
    required_vars = exprs(PREFIX, GRPNAME, SRCVAR, TERMNAME, TERMID)
  )

  # check duplicate rows
  signal_duplicate_records(queries, by_vars = exprs(!!!syms(colnames(queries))))

  # check illegal prefix category
  is_good_prefix <- grepl("^[a-zA-Z]{2,3}", queries$PREFIX)
  if (!all(is_good_prefix)) {
    abort(
      paste0(
        "`PREFIX` in `", queries_name,
        "` must start with 2-3 letters.. Problem with ",
        enumerate(unique(queries$PREFIX[!is_good_prefix])),
        "."
      )
    )
  }

  # check illegal prefix number
  query_num <- sub("[[:alpha:]]+", "", queries$PREFIX)
  is_bad_num <- nchar(query_num) != 2 | is.na(as.numeric(query_num))
  if (any(is_bad_num)) {
    abort(
      paste0(
        "`PREFIX` in `", queries_name,
        "` must end with 2-digit numbers. Issue with ",
        enumerate(unique(queries$PREFIX[is_bad_num])),
        "."
      )
    )
  }

  # check illegal query name
  if (any(queries$GRPNAME == "") || any(is.na(queries$GRPNAME))) {
    abort(paste0(
      "`GRPNAME` in `", queries_name,
      "` cannot be empty string or NA."
    ))
  }

  # check query id is numeric
  if ("GRPID" %in% names(queries) && !is.numeric(queries$GRPID)) {
    abort(paste0(
      "`GRPID` in `", queries_name,
      "` should be numeric."
    ))
  }

  # check illegal query scope
  if ("SCOPE" %in% names(queries) &&
    any(unique(queries$SCOPE) %notin% c("BROAD", "NARROW", "", NA_character_))) {
    abort(paste0(
      "`SCOPE` in `", queries_name,
      "` can only be 'BROAD', 'NARROW' or `NA`."
    ))
  }

  # check illegal query scope number
  if ("SCOPEN" %in% names(queries)) {
    is_bad_scope_num <- queries$SCOPEN %notin% c(1, 2, NA_integer_)
    if (any(is_bad_scope_num)) {
      abort(
        paste0(
          "`SCOPEN` in `", queries_name,
          "` must be one of 1, 2, or NA. Issue with ",
          enumerate(unique(queries$SCOPEN[is_bad_scope_num])),
          "."
        )
      )
    }
  }

  # check illegal term name
  if (any(is.na(queries$TERMNAME) & is.na(queries$TERMID)) ||
    any(queries$TERMNAME == "" & is.na(queries$TERMID))) {
    abort(paste0(
      "Either `TERMNAME` or `TERMID` need to be specified",
      " in `", queries_name, "`. ",
      "They both cannot be NA or empty."
    ))
  }

  # each PREFIX must have unique GRPNAME, GRPID if the columns exist
  count_unique <- queries %>%
    group_by(PREFIX) %>%
    summarise(
      n_qnam = length(unique(GRPNAME)),
      n_qid = ifelse("GRPID" %in% names(queries),
        length(unique(GRPID)), 0
      )
    ) %>%
    ungroup()

  if (any(count_unique$n_qnam > 1)) {
    idx <- which(count_unique$n_qnam > 1)
    abort(paste0(
      "In `", queries_name, "`, `GRPNAME` of '",
      paste(count_unique$PREFIX[idx], collapse = ", "),
      "' is not unique."
    ))
  }

  if (any(count_unique$n_qid > 1)) {
    idx <- which(count_unique$n_qid > 1)
    abort(paste0(
      "In `", queries_name, "`, `GRPID` of '",
      paste(count_unique$PREFIX[idx], collapse = ", "),
      "' is not unique."
    ))
  }

  # check SCOPE and SCOPEN are one to one if available
  if ("SCOPE" %in% names(queries) && "SCOPEN" %in% names(queries)) {
    assert_one_to_one(queries, exprs(SCOPE), exprs(SCOPEN))
  }
}
