#' Derive Variables by Transposing and Merging a Second Dataset
#'
#' Adds variables from a vertical dataset after transposing it into a wide one.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_merge Dataset to transpose and merge
#'
#'   The variables specified by the `by_vars`, `id_vars`, `key_var` and
#'   `value_var` arguments are expected. The variables `by_vars`, `id_vars`,
#'   `key_var` have to be a unique key.
#'
#' @param by_vars Grouping variables
#'
#'  Keys used to merge `dataset_merge` with `dataset`.
#'
#' @param id_vars ID variables
#'
#'  Variables (excluding by_vars) that uniquely identify each observation in `dataset_merge`.
#'
#' `r roxygen_param_by_vars()`
#'
#' @param key_var The variable of `dataset_merge` containing the names of the
#'   transposed variables
#'
#' @param value_var The variable of `dataset_merge` containing the values of the
#'   transposed variables
#'
#' @param filter Expression used to restrict the records of `dataset_merge` prior to transposing
#'
#' @param relationship Expected merge-relationship between the `by_vars`
#'   variable(s) in `dataset` and `dataset_merge` (after transposition)
#'
#'   This argument is passed to the `dplyr::left_join()` function. See
#'   <https://dplyr.tidyverse.org/reference/mutate-joins.html#arguments> for
#'   more details.
#'
#'   *Permitted Values*: `"one-to-one"`, `"one-to-many"`, `"many-to-one"`,
#'   `"many-to-many"`, `NULL`
#'
#' @details
#' After filtering `dataset_merge` based upon the condition provided in `filter`, this
#' dataset is transposed and subsequently merged onto `dataset` using `by_vars` as
#' keys.
#'
#' @return The input dataset with transposed variables from `dataset_merge` added
#'
#' @seealso [derive_vars_atc()]
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' cm <- tribble(
#'   ~USUBJID,       ~CMGRPID, ~CMREFID,  ~CMDECOD,
#'   "BP40257-1001", "14",     "1192056", "PARACETAMOL",
#'   "BP40257-1001", "18",     "2007001", "SOLUMEDROL",
#'   "BP40257-1002", "19",     "2791596", "SPIRONOLACTONE"
#' )
#' facm <- tribble(
#'   ~USUBJID,       ~FAGRPID, ~FAREFID,  ~FATESTCD,  ~FASTRESC,
#'   "BP40257-1001", "1",      "1192056", "CMATC1CD", "N",
#'   "BP40257-1001", "1",      "1192056", "CMATC2CD", "N02",
#'   "BP40257-1001", "1",      "1192056", "CMATC3CD", "N02B",
#'   "BP40257-1001", "1",      "1192056", "CMATC4CD", "N02BE",
#'   "BP40257-1001", "1",      "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "1",      "2007001", "CMATC2CD", "D10",
#'   "BP40257-1001", "1",      "2007001", "CMATC3CD", "D10A",
#'   "BP40257-1001", "1",      "2007001", "CMATC4CD", "D10AA",
#'   "BP40257-1001", "2",      "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "2",      "2007001", "CMATC2CD", "D07",
#'   "BP40257-1001", "2",      "2007001", "CMATC3CD", "D07A",
#'   "BP40257-1001", "2",      "2007001", "CMATC4CD", "D07AA",
#'   "BP40257-1001", "3",      "2007001", "CMATC1CD", "H",
#'   "BP40257-1001", "3",      "2007001", "CMATC2CD", "H02",
#'   "BP40257-1001", "3",      "2007001", "CMATC3CD", "H02A",
#'   "BP40257-1001", "3",      "2007001", "CMATC4CD", "H02AB",
#'   "BP40257-1002", "1",      "2791596", "CMATC1CD", "C",
#'   "BP40257-1002", "1",      "2791596", "CMATC2CD", "C03",
#'   "BP40257-1002", "1",      "2791596", "CMATC3CD", "C03D",
#'   "BP40257-1002", "1",      "2791596", "CMATC4CD", "C03DA"
#' )
#'
#' cm %>%
#'   derive_vars_transposed(
#'     facm,
#'     by_vars = exprs(USUBJID, CMREFID = FAREFID),
#'     id_vars = exprs(FAGRPID),
#'     key_var = FATESTCD,
#'     value_var = FASTRESC
#'   ) %>%
#'   select(USUBJID, CMDECOD, starts_with("CMATC"))
derive_vars_transposed <- function(dataset,
                                   dataset_merge,
                                   by_vars,
                                   id_vars = NULL,
                                   key_var,
                                   value_var,
                                   filter = NULL,
                                   relationship = NULL) {
  key_var <- assert_symbol(enexpr(key_var))
  value_var <- assert_symbol(enexpr(value_var))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_vars(by_vars)
  assert_vars(id_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(dataset_merge, required_vars = expr_c(by_vars, key_var, value_var))
  relationship <- assert_character_scalar(
    relationship,
    values = c("one-to-one", "one-to-many", "many-to-one", "many-to-many"),
    case_sensitive = TRUE,
    optional = TRUE
  )

  dataset_merge <- filter_if(dataset_merge, filter)

  # check for duplicates in dataset_merge as these will create list columns,
  # which is not acceptable for ADaM datasets
  signal_duplicate_records(
    dataset_merge,
    by_vars = c(by_vars, id_vars, exprs(!!key_var)),
    msg = c(
      paste(
        "Dataset {.arg dataset_merge} contains duplicate records with respect to",
        "{.var {by_vars}}"
      ),
      "Please check data and {.arg by_vars}, {.arg id_vars}, and {.arg key_var} arguments."
    ),
    class = "merge_duplicates"
  )

  dataset_transposed <- dataset_merge %>%
    pivot_wider(
      names_from = !!key_var,
      values_from = !!value_var,
      id_cols = c(as.character(by_vars), as.character(id_vars))
    )

  tryCatch(
    left_join(
      dataset,
      dataset_transposed,
      by = vars2chr(by_vars),
      relationship = relationship
    ),
    "dplyr_error_join_relationship_one_to_one" = function(cnd) {
      cli_abort(
        message = c(
          str_replace(
            str_replace(
              cnd$message, "`x`", "`dataset`"
            ), "`y`", "the transposed `dataset_merge`"
          ),
          i = str_replace(
            str_replace(
              cnd$body, "`x`", "`dataset`"
            ), "`y`", "the transposed `dataset_merge`"
          )
        ),
        call = parent.frame(n = 4)
      )
    },
    "dplyr_error_join_relationship_many_to_one" = function(cnd) {
      cli_abort(
        message = c(
          str_replace(
            str_replace(
              cnd$message, "`x`", "`dataset`"
            ), "`y`", "the transposed `dataset_merge`"
          ),
          i = str_replace(
            str_replace(
              cnd$body, "`x`", "`dataset`"
            ), "`y`", "the transposed `dataset_merge`"
          )
        ),
        call = parent.frame(n = 4)
      )
    }
  )
}
