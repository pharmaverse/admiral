filter_extreme_multisource <- function(sources,
                                       by_vars,
                                       order,
                                       mode){
  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (!is.null(sources[[i]]$filter)) {
      add_data[[i]] <- sources[[i]]$dataset %>%
        filter(!!(sources[[i]]$filter))
    }
    else {
      add_data[[i]] <- sources[[i]]$dataset
    }
    if (!is.null(sources[[i]]$order)){
      add_data[[i]] <- filter_extreme(add_data[[i]],
                                      order = sources[[i]]$order,
                                      by_vars = by_vars,
                                      mode = sources[[i]]$mode)
    }
    add_data[[i]] <- transmute(add_data[[i]], temp_source_nr = i, !!!by_vars, !!!sources[[i]]$vars)
  }
  bind_rows(add_data) %>%
    filter_extreme(by_vars = by_vars,
                   order = order,
                   mode = mode) %>%
    select(-starts_with("temp_"))
}
