#' Get Query Variables
#'
#' @description Create a table for the input dataset which binds the necessary
#' rows for a `derive_vars_query()` call with the relevant `SRCVAR`, `TERM_NAME_ID`
#' and a temporary index if it is necessary
#'
#' **Note:** This function is the first step performed in `derive_vars_query()`
#' requested by some users to be present independently from it.
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
#'   `GRPNAME` if the value of `TERMCHAR` or `TERMNUM` in `dataset_queries` matches
#'   the value of the respective SRCVAR in `dataset`.
#'   Note that `TERMCHAR` in `dataset_queries` dataset may be NA only when `TERMNUM`
#'   is non-NA and vice versa. The matching is case insensitive.
#'   The "CD", "SC", and "SCN" variables are derived accordingly based on
#'   `GRPID`, `SCOPE`, and `SCOPEN` respectively,
#'   whenever not missing.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' @param dataset_queries A dataset containing required columns `PREFIX`,
#' `GRPNAME`, `SRCVAR`, `TERMCHAR` and/or `TERMNUM`, and optional columns
#' `GRPID`, `SCOPE`, `SCOPEN`.
#'
#' `create_query_data()` can be used to create the dataset.
#'
#'
#' @return The processed query dataset with `SRCVAR` and `TERM_NAME_ID` so that
#' that can be merged to the input dataset to execute the derivations outlined by `dataset_queries`.
#'
#' @family utils_help
#' @keywords utils_help
#'
#' @seealso [create_query_data()]
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
#' get_vars_query(adae, queries)
get_vars_query <- function(dataset, dataset_queries) { # nolint: cyclocomp_linter
  source_vars <- unique(dataset_queries$SRCVAR)
  assert_data_frame(dataset,
    required_vars = chr2vars(source_vars),
    optional = FALSE
  )

  # check optionality of TERMNUM or TERMCHAR based on SRCVAR type
  srcvar_types <- unique(vapply(dataset[source_vars], typeof, character(1)))
  if (!all(srcvar_types %in% c("character", "integer", "double"))) {
    idx <- source_vars[!vapply(dataset[source_vars], typeof, character(1)) %in% c("character", "integer", "double")] # nolint
    dat_incorrect_type <- dataset[idx]
    cli_abort(c(
      "The source variables (values of {.var SRCVAR}) must be numeric or character.",
      i = paste0(
        colnames(dat_incorrect_type),
        " is of type ",
        vapply(dat_incorrect_type, typeof, character(1)),
        collapse = ", "
      )
    ))
  }

  termvars <- exprs(character = TERMCHAR, integer = TERMNUM, double = TERMNUM)
  expected_termvars <- unique(termvars[srcvar_types])
  assert_data_frame(dataset_queries, required_vars = c(exprs(PREFIX, GRPNAME, SRCVAR), expected_termvars)) # nolint
  if (length(expected_termvars) > 1) {
    # check illegal term name
    if (any(is.na(dataset_queries$TERMCHAR) & is.na(dataset_queries$TERMNUM)) ||
      any(dataset_queries$TERMCHAR == "" & is.na(dataset_queries$TERMNUM))) {
      cli_abort(paste0(
        "Either {.var TERMCHAR} or {.var TERMNUM} need to be specified",
        " in {.arg dataset_queries}. ",
        "They both cannot be NA or empty."
      ))
    }
  }
  assert_valid_queries(dataset_queries, queries_name = deparse(substitute(dataset_queries)))

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
    # numeric -> TERMNUM, character -> TERMCHAR, otherwise -> error
    mutate(
      tmp_col_type = vapply(dataset[SRCVAR], typeof, character(1)),
      TERM_NAME_ID = ifelse(
        tmp_col_type == "character",
        toupper(TERMCHAR),
        as.character(TERMNUM)
      )
    )

  # prepare input dataset for joining
  static_cols <- setdiff(names(dataset), chr2vars(source_vars))
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
  joined %>%
    inner_join(queries_wide, by = c("SRCVAR", "TERM_NAME_ID")) %>%
    select(!!!syms(c(static_cols, new_col_names))) %>%
    group_by_at(static_cols) %>%
    summarise_all(~ first(na.omit(.))) %>%
    ungroup()
}

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
#'   `GRPNAME` if the value of `TERMCHAR` or `TERMNUM` in `dataset_queries` matches
#'   the value of the respective SRCVAR in `dataset`.
#'   Note that `TERMCHAR` in `dataset_queries` dataset may be NA only when `TERMNUM`
#'   is non-NA and vice versa. The matching is case insensitive.
#'   The "CD", "SC", and "SCN" variables are derived accordingly based on
#'   `GRPID`, `SCOPE`, and `SCOPEN` respectively,
#'   whenever not missing.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' @param dataset_queries A dataset containing required columns `PREFIX`,
#' `GRPNAME`, `SRCVAR`, `TERMCHAR` and/or `TERMNUM`, and optional columns
#' `GRPID`, `SCOPE`, `SCOPEN`.
#'
#' `create_query_data()` can be used to create the dataset.
#'
#'
#' @return The input dataset with query variables derived.
#'
#' @family der_occds
#' @keywords der_occds
#'
#' @seealso [create_query_data()]
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
derive_vars_query <- function(dataset, dataset_queries) { # nolint: cyclocomp_linter
  # join restructured queries to input dataset
  assert_valid_queries(dataset_queries, queries_name = deparse(substitute(dataset_queries)))
  dataset_queries <- convert_blanks_to_na(dataset_queries)
  source_vars <- unique(dataset_queries$SRCVAR)
  static_cols <- setdiff(names(dataset), chr2vars(source_vars))
  no_key <- dataset %>%
    select(all_of(static_cols)) %>%
    distinct()
  if (nrow(no_key) != nrow(dataset)) {
    dataset$temp_key <- seq_len(nrow(dataset))
    static_cols <- c(static_cols, "temp_key")
  }
  tryCatch(
    expr = {
      joined <- get_vars_query(dataset, dataset_queries)
    },
    error = function(e) {
      stop("Error in derive_vars_query call of get_vars_query: ", e)
    }
  )
  # join queries to input dataset, remove temp col(s)
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
#' - `TERMCHAR`, character, could be NA only at those observations
#' where `TERMNUM` is non-NA
#' - `TERMNUM`, integer, could be NA only at those observations
#' where `TERMCHAR` is non-NA
#'
#' @param queries A data.frame.
#'
#' @param queries_name Name of the queries dataset, a string.
#'
#' @return The function throws an error if any of the requirements not met.
#'
#' @examples
#' data("queries")
#' assert_valid_queries(queries, "queries")
#' @noRd
assert_valid_queries <- function(queries, queries_name) {
  # check duplicate rows
  signal_duplicate_records(queries, by_vars = exprs(!!!syms(colnames(queries))))

  # check illegal prefix category
  is_good_prefix <- grepl("^[a-zA-Z]{2,3}", queries$PREFIX)
  if (!all(is_good_prefix)) {
    cli_abort(
      paste0(
        "{.var PREFIX} in {.arg {queries_name}}",
        " must start with 2-3 letters. Problem with ",
        "{.val {unique(queries$PREFIX[!is_good_prefix])}}."
      )
    )
  }

  # check illegal prefix number
  query_num <- sub("[[:alpha:]]+", "", queries$PREFIX)
  is_bad_num <- nchar(query_num) != 2 | is.na(as.numeric(query_num))
  if (any(is_bad_num)) {
    cli_abort(
      paste0(
        "{.var PREFIX} in {.arg {queries_name}}",
        " must end with 2-digit numbers. Issue with ",
        "{.val {unique(queries$PREFIX[is_bad_num])}}."
      )
    )
  }

  # check illegal query name
  if (any(queries$GRPNAME == "") || any(is.na(queries$GRPNAME))) {
    cli_abort(
      "{.var GRPNAME} in {.arg {queries_name}} cannot be empty string or NA."
    )
  }

  # check query id is numeric
  if ("GRPID" %in% names(queries) && !is.numeric(queries$GRPID)) {
    cli_abort(
      "{.var GRPID} in {.arg {queries_name}} must be numeric."
    )
  }

  # check illegal query scope
  if ("SCOPE" %in% names(queries) &&
    any(unique(queries$SCOPE) %notin% c("BROAD", "NARROW", "", NA_character_))) {
    cli_abort(
      "{.var SCOPE} in {.arg {queries_name}} can only be 'BROAD', 'NARROW' or `NA`."
    )
  }

  # check illegal query scope number
  if ("SCOPEN" %in% names(queries)) {
    is_bad_scope_num <- queries$SCOPEN %notin% c(1, 2, NA_integer_)
    if (any(is_bad_scope_num)) {
      cli_abort(
        paste0(
          "{.var SCOPEN} in {.arg {queries_name}}",
          " must be one of 1, 2, or NA. Issue with ",
          "{.val {unique(queries$SCOPEN[is_bad_scope_num])}}."
        )
      )
    }
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
    cli_abort(paste0(
      "In {.arg {queries_name}} {.var GRPNAME} of ",
      "{.val {count_unique$PREFIX[idx]}} is not unique."
    ))
  }

  if (any(count_unique$n_qid > 1)) {
    idx <- which(count_unique$n_qid > 1)
    cli_abort(paste0(
      "In {.arg {queries_name}} {.var GRPID} of ",
      "{.val {count_unique$PREFIX[idx]}} is not unique."
    ))
  }

  # check SCOPE and SCOPEN are one to one if available
  if ("SCOPE" %in% names(queries) && "SCOPEN" %in% names(queries)) {
    assert_one_to_one(queries, exprs(SCOPE), exprs(SCOPEN))
  }
}
