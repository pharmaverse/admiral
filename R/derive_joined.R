derive_vars_joined <- function(dataset,
                               dataset_add,
                               by_vars = NULL,
                               order = NULL,
                               new_vars = NULL,
                               join_vars,
                               filter_add = NULL,
                               filter_join,
                               mode = NULL) {
  assert_vars(by_vars, optional = TRUE)
  assert_order_vars(order, optional = TRUE)
  assert_vars(new_vars, optional = TRUE)
  assert_vars(join_vars)
  assert_data_frame(dataset, required_vars = by_vars)
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, join_vars, extract_vars(order), new_vars))
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enquo(filter_join))

  if (is.null(new_vars)) {
    new_vars = quos(!!!syms(colnames(dataset_add)))
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

  data_joined <- left_join(
    data,
    filter_if(dataset_add, filter_add),
    by = vars2chr(by_vars),
    suffix = c("", ".join")
  )

  data_return <- filter(data_joined, !!filter_join)

  if (!is.null(order)) {
    data_return <- filter_extreme(
      data_return,
      by_vars = quo_c(by_vars, tmp_obs_nr),
      order = order,
      mode = mode
    )
  }
  data %>%
    derive_vars_merged(
      dataset_add = select(data_return, !!tmp_obs_nr, !!!new_vars),
      by_vars = vars(!!tmp_obs_nr)
    ) %>%
    remove_tmp_vars()
}
