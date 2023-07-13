# derive_var_basetype ----
## Test 1: deprecation error if function is called ----
test_that("derive_var_basetype Test 1: deprecation error if function is called", {
  input <- tibble::tribble(
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
  expect_output <- tibble::tribble(
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
  expect_error(
    derive_var_basetype(
      dataset = input,
      basetypes = rlang::exprs(
        "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
        "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
        "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
      )
    ),
    class = "lifecycle_error_deprecated"
  )
})
