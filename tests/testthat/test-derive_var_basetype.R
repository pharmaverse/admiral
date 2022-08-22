test_that("records are duplicated across different `BASETYPE` values", {
  library(tibble)
  input <- tribble(
    ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
    "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
    "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
    "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1,
    "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4,
    "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9,
    "P02",    "RUN-IN",       "PARAM01",     1,  12.1,
    "P02",    "DOUBLE-BLIND", "PARAM01",     2,  10.2,
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.8,
    "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4,
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8
  )
  expect_output <- tribble(
    ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL, ~BASETYPE,
    "P01",    "RUN-IN",       "PARAM01",     1,  10.0, "RUN-IN",
    "P01",    "RUN-IN",       "PARAM01",     2,   9.8, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1, "RUN-IN",
    "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4, "RUN-IN",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2, "DOUBLE-BLIND",
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4, "OPEN-LABEL",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9, "OPEN-LABEL",
    "P02",    "RUN-IN",       "PARAM01",     1,  12.1, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     2,  10.2, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.8, "RUN-IN",
    "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4, "RUN-IN",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     2,  10.2, "DOUBLE-BLIND",
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.8, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4, "OPEN-LABEL",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8, "OPEN-LABEL",
  )
  actual_output <- derive_var_basetype(
    dataset = input,
    basetypes = rlang::exprs(
      "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
      "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
      "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
    )
  )

  expect_dfs_equal(actual_output, expect_output, keys = c("USUBJID", "BASETYPE", "PARAMCD", "ASEQ"))
})

test_that("records that do not match any condition are kept", {
  library(tibble)
  input <- tribble(
    ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
    "P01",    "SCREENING",    "PARAM01",     1,  10.2,
    "P01",    "RUN-IN",       "PARAM01",     2,  10.0,
    "P01",    "RUN-IN",       "PARAM01",     3,   9.8,
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,   9.2,
    "P01",    "DOUBLE-BLIND", "PARAM01",     5,  10.1,
    "P01",    "OPEN-LABEL",   "PARAM01",     6,  10.4,
    "P01",    "OPEN-LABEL",   "PARAM01",     7,   9.9,
    "P02",    "SCREENING",    "PARAM01",     1,  12.2,
    "P02",    "RUN-IN",       "PARAM01",     2,  12.1,
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.2,
    "P02",    "DOUBLE-BLIND", "PARAM01",     4,  10.8,
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  11.4,
    "P02",    "OPEN-LABEL",   "PARAM01",     6,  10.8
  )
  expect_output <- tribble(
    ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL, ~BASETYPE,
    "P01",    "SCREENING",    "PARAM01",     1,  10.2, NA,
    "P01",    "RUN-IN",       "PARAM01",     2,  10.0, "RUN-IN",
    "P01",    "RUN-IN",       "PARAM01",     3,   9.8, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,   9.2, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     5,  10.1, "RUN-IN",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,  10.4, "RUN-IN",
    "P01",    "OPEN-LABEL",   "PARAM01",     7,   9.9, "RUN-IN",
    "P01",    "DOUBLE-BLIND", "PARAM01",     4,   9.2, "DOUBLE-BLIND",
    "P01",    "DOUBLE-BLIND", "PARAM01",     5,  10.1, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,  10.4, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     7,   9.9, "DOUBLE-BLIND",
    "P01",    "OPEN-LABEL",   "PARAM01",     6,  10.4, "OPEN-LABEL",
    "P01",    "OPEN-LABEL",   "PARAM01",     7,   9.9, "OPEN-LABEL",
    "P02",    "SCREENING",    "PARAM01",     1,  12.2, NA,
    "P02",    "RUN-IN",       "PARAM01",     2,  12.1, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.2, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     4,  10.8, "RUN-IN",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  11.4, "RUN-IN",
    "P02",    "OPEN-LABEL",   "PARAM01",     6,  10.8, "RUN-IN",
    "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.2, "DOUBLE-BLIND",
    "P02",    "DOUBLE-BLIND", "PARAM01",     4,  10.8, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  11.4, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     6,  10.8, "DOUBLE-BLIND",
    "P02",    "OPEN-LABEL",   "PARAM01",     5,  11.4, "OPEN-LABEL",
    "P02",    "OPEN-LABEL",   "PARAM01",     6,  10.8, "OPEN-LABEL",
  )
  actual_output <- derive_var_basetype(
    dataset = input,
    basetypes = rlang::exprs(
      "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
      "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
      "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
    )
  )

  expect_dfs_equal(actual_output, expect_output, keys = c("USUBJID", "BASETYPE", "PARAMCD", "ASEQ"))
})
