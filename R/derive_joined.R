derive_vars_joined <- function(dataset,
                               dataset_add,
                               by_vars = NULL,
                               order = NULL,
                               new_vars = NULL,
                               join_vars,
                               filter_add = NULL,
                               filter_join,
                               mode = NULL,
                               check_type = "warning") {
  assert_vars(by_vars, optional = TRUE)
  assert_order_vars(order, optional = TRUE)
  assert_vars(new_vars, optional = TRUE)
  assert_vars(join_vars)
  assert_data_frame(dataset, required_vars = by_vars)
#  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, join_vars, extract_vars(order), new_vars))
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enquo(filter_join))

  if (is.null(new_vars)) {
    new_vars = chr2vars(colnames(dataset_add))
  }

  # number observations of the input dataset to get a unique key
  # (by_vars and tmp_obs_nr)
  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr,
      by_vars = by_vars,
      check_type = "none"
    )

  data_right <- filter_if(dataset_add, filter_add) %>%
    select(!!!by_vars, !!!join_vars, !!!unname(new_vars))
  data_joined <- left_join(
    data,
    data_right,
    by = vars2chr(by_vars),
    suffix = c("", ".join")
  )

  data_return <- filter(data_joined, !!filter_join)

  common_vars <- chr2vars(intersect(colnames(data), colnames(data_right)))
  if (!is.null(order)) {
    data_return <- filter_extreme(
      data_return,
      by_vars = quo_c(by_vars, quo(!!tmp_obs_nr)),
      order = add_suffix_to_vars(order, vars = common_vars, suffix = ".join"),
      mode = mode,
      check_type = check_type
    )
  }
  data %>%
    derive_vars_merged(
      dataset_add = select(
        data_return,
        !!tmp_obs_nr,
        !!!add_suffix_to_vars(new_vars, vars = common_vars, suffix = ".join")
      ),
      by_vars = vars(!!tmp_obs_nr)
    ) %>%
    remove_tmp_vars()
}


replace_symbol_in_quo <- function(quosure,
                                  target,
                                  replace) {
  assert_expr(quosure)
  target <- quo_get_expr(assert_symbol(enquo(target)))
  replace <- quo_get_expr(assert_symbol(enquo(replace)))
  expr <- quo_get_expr(quosure)
  if (is.symbol(expr)) {
    if (expr == target) {
      expr = replace
    }
  } else {
    for (i in seq_along(quosure)) {
      if (expr[[i]] == target) {
        expr[[i]] <- replace
      }
    }
  }
  rlang::quo_set_expr(quosure, expr)
}

add_suffix_to_vars <- function(order,
                               vars,
                               suffix
) {
  assert_order_vars(order)
  assert_vars(vars)
  assert_character_scalar(suffix)
  for (i in seq_along(vars)) {
    order <- lapply(
      order,
      replace_symbol_in_quo,
      target = !!quo_get_expr(vars[[i]]),
      replace = !!sym(paste0(as_label(
             quo_get_expr(vars[[i]])
           ), suffix)))
  }
  class(order) <- c("quosures", "list")
  order
}
