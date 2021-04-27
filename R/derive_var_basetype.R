derive_var_basetype <- function(dataset, basetypes) {
  assert_that(
    is.data.frame(dataset),
    is_named_exprs(basetypes)
  )
  assert_has_variables(
    dataset,
    unique(map_chr(basetypes, all.vars))
  )

  subsets <- map2(names(basetypes), basetypes, function(label, condition) {
    dataset %>%
      filter(!!condition) %>%
      mutate(BASETYPE = label)
  })

  bind_rows(subsets)
}
