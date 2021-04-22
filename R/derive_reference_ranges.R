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

  has_a1lo <- "A1LO" %in% colnames(meta_ref_ranges)
  has_a1hi <- "A1HI" %in% colnames(meta_ref_ranges)
  if (!has_a1lo) meta_ref_ranges$A1LO <- NA
  if (!has_a1hi) meta_ref_ranges$A1HI <- NA

  result <- dataset %>%
    left_join(meta_ref_ranges, by = by_var) %>%
    mutate(
      ANRIND = case_when(
        AVAL >= ANRLO & AVAL <= ANRHI ~ "NORMAL",
        AVAL < ANRLO & (is.na(A1LO) | AVAL >= A1LO) ~ "LOW",
        AVAL > ANRHI & (is.na(A1HI) | AVAL <= A1HI) ~ "HIGH",
        AVAL < A1LO ~ "LOW LOW",
        AVAL > A1HI ~ "HIGH HIGH",
        TRUE ~ ""
      )
    )

  if (!has_a1lo) result$A1LO <- NULL
  if (!has_a1hi) result$A1HI <- NULL

  result
}
