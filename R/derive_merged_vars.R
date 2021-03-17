
derive_merged_vars <- function(dataset,
                               dataset_add,
                               filter_add,
                               new_vars,
                               by_vars = exprs(USUBJID),
                               filter_first_order){
  if (!missing(filter_add)){
     add <- dataset_add %>%
       filter(!!!filter_add)
  }
  else{
    add <- dataset_add
  }
  if (!missing(filter_first_order)){
    add <- add %>%
      filter_first(order = filter_first_order,
                   by_vars = by_vars)
  }
  add <- add %>% transmute(!!!new_vars)
  result <- left_join(dataset, add, by = map_chr(by_vars, as_string))
  return(result)
}
