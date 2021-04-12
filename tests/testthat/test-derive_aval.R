context("test-derive_aval")

test_that("`--STRES`, `--STRESN` and `--STRESU` are mapped to `AVALC`, `AVAL` and `AVALU`", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~VSSTRES, ~VSSTRESN, ~VSSTRESU,
    "P01",    "HEIGHT", "180.3",  180.3,     "cm",
    "P02",    "HEIGHT", "159.9",  159.9,     "cm",
    "P03",    "HEIGHT", "167.0",  167,       "cm"
  )
  new_cols <- tibble::tribble(
    ~AVALC,  ~AVAL, ~AVALU,
    "180.3", 180.3, "cm",
    "159.9", 159.9, "cm",
    "167.0", 167,   "cm"
  )
  expected_output <- dplyr::bind_cols(input, new_cols)

  expect_dfs_equal(expected_output, derive_aval(input), keys = c("USUBJID", "PARAMCD"))
})
