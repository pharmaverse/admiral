library(admiral.test)

test_that("a warning is issued when specifying `derive_extreme_flag(flag_filter = )`", {
  data(advs)

  expect_warning(
    derive_extreme_flag(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAMCD),
      order = vars(ADT),
      new_var = ABLFL,
      mode = "last",
      flag_filter = AVISIT == "BASELINE"
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `dthcaus_source(date_var = )", {
  expect_warning(
    dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date_var = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `dthcaus_source(traceabilty_vars = )", {
  expect_warning(
    dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD,
      traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `lstalvdt_source(date_var = )", {
  expect_warning(
    lstalvdt_source(dataset_name = "adsl", date_var = TRTEDT),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_suppqual_vars()", {
  data(ae)
  data(suppae)

  expect_warning(
    derive_suppqual_vars(ae[1:100, ], suppae[1:100, ]),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_query_vars()", {
  data(queries)
  adae <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~AELLTCD,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
      3, "Alanine aminotransferase abnormal", NA_character_, NA_integer_,
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
      5, "Basedow's disease", NA_character_, 1L,
    "03", "2020-06-07 23:59:59", "SOME TERM",
      2, "Some query", "Some term", NA_integer_,
    "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
      7, "Alveolar proteinosis", NA_character_,  NA_integer_
  )

  expect_warning(
    derive_query_vars(adae, queries),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_duration()", {
  adsl <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
  )

  expect_warning(
    derive_duration(
      adsl,
      new_var = AAGE,
      new_var_unit = AAGEU,
      start_date = BRTHDT,
      end_date = RANDDT,
      out_unit = "years",
      add_one = FALSE,
      trunc_out = TRUE
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_aage()", {
  adsl <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
  )

  expect_warning(
    derive_aage(adsl),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `derive_var_ontrtfl(date = )", {
  data(advs)

  expect_warning(
    derive_var_ontrtfl(
      advs[1:100, ],
      date = ADT,
      ref_start_date = TRTSDT,
      ref_end_date = TRTEDT
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `derive_summary_records(filter_rows = )", {
  data(advs)

  expect_warning(
    derive_summary_records(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAM, AVISIT),
      analysis_var = AVAL,
      summary_fun = function(x) mean(x, na.rm = TRUE),
      filter_rows = dplyr::n() > 2,
      set_values_to = vars(DTYPE = "AVERAGE")
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("an error is thrown when specifying `derive_summary_records(fns = )", {
  data(advs)

  expect_error(
    derive_summary_records(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAM, AVISIT),
      fns = AVAL ~ mean(., na.rm = TRUE),
      filter = dplyr::n() > 2,
      set_values_to = vars(DTYPE = "AVERAGE")
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `dthcaus_source(dataset = )", {
  expect_warning(
    dthcaus_source(
      dataset = ae,
      filter = AEOUT == "FATAL",
      date = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `lstalvdt_source(dataset = )", {
  expect_warning(
    lstalvdt_source(
      dataset = lb,
      date = LBDTC,
      filter = nchar(LBDTC) >= 10
    ),
    "deprecated",
    fixed = TRUE
  )
})
