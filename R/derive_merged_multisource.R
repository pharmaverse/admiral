derive_merged_multisource <- function(dataset,
                                      sources,
                                      by_vars,
                                      order,
                                      mode){
  all_data <- filter_extreme_multisource(sources = sources,
                                         by_vars = by_vars,
                                         order = order,
                                         mode = mode)
  left_join(dataset, all_data, by = map_chr(by_vars, as_string)) %>%
    select(-starts_with("temp_"))
}
