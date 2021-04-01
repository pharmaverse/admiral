derive_obs_number <- function(dataset,
                              new_var = temp_obs_nr ,
                              order,
                              by_vars,
                              check_type = 'warning'){
  arg_match(check_type, c('none', 'warning', 'error'))
  data <- dataset

  if (!missing(by_vars) | !missing(order)) {
    # group and sort input dataset
    if (!missing(by_vars)) {
      assert_has_variables(dataset, by_vars)

      data <- dataset %>% group_by(!!!by_vars) %>%
        arrange(!!!order, .by_group = TRUE)

      if (check_type != 'none') {
        # check for unique records
        assert_has_unique_records(
          data,
          by_vars = by_vars ,
          order = order ,
          message_type = check_type
        )
      }
    }
    else{
      data <- dataset %>%
        arrange(!!!order)

      if (check_type != 'none') {
        # check for unique records
        assert_has_unique_records(data,
                                  order = order,
                                  message_type = 'warning')
      }
    }
  }

  data %>%
    mutate(!!enquo(new_var) := row_number()) %>%
    ungroup()
}
