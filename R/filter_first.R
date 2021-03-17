filter_first <- function(dataset,
                         order,
                         by_vars){
  dataset %>% group_by(!!!by_vars) %>%
    arrange(!!!order, .by_group = TRUE) %>%
    slice(1)
}
