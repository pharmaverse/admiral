derive_var_basetype <- function(dataset, basetypes) {
  subsets <- map2(names(basetypes), basetypes, function(label, condition) {
    dataset %>%
      filter(!!condition) %>%
      mutate(BASETYPE = label)
  })

  bind_rows(subsets)
}
