derive_reference_ranges <- function(dataset,
                                    meta_ref_ranges,
                                    by_var = "PARAMCD") {
  assert_that(
    is.data.frame(dataset),
    is.data.frame(meta_ref_ranges),
    is.character(by_var)
  )
  assert_has_variables(dataset, by_var)
  assert_has_variables(meta_ref_ranges, c(by_var, "ANRLO", "ANRHI"))

  warn_if_ref_ranges_missing(dataset, meta_ref_ranges, by_var)

  left_join(dataset, meta_ref_ranges, by = by_var)
}
