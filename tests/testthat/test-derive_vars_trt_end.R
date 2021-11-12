context("test-derive_vars_trt_end")

#test 1
test_that("Test1::Adds treatment end date - no imputations
          choosing the LAST ex date with filter of EXDOSE>0 applied", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06-19", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDTM,
              "TEST01", "PAT01", "2014-01-17T00:00:00",
              "TEST01", "PAT02", "2012-06-19T00:00:00"
            ) %>% mutate(
              TRTEDTM = as_iso_dtm(TRTEDTM),
            )
            # nolint end

            actual_output <- derive_vars_trt_end(input,
                                                   input_ex,
                                                   new_var=TRTEDTM,
                                                   dtc=EXSTDTC,
                                                   date_imputation = "first",
                                                   time_imputation="first",
                                                   min_dates = NULL,
                                                   max_dates = NULL,
                                                   order = vars(EXSTDTC,EXSEQ),
                                                   subject_keys=vars(STUDYID,USUBJID),
                                                   flag_imputation = "AUTO",
                                                   filter_ex = EXDOSE>0 ,
                                                   ord_filter="LAST")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

#test 2
test_that("Test2::Adds treatment start date - no imputations
          choosing the LAST ex date with no filter applied
          -expect a warning thrown for no filter applied", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06-19", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDTM,
              "TEST01", "PAT01", "2014-06-19T00:00:00",
              "TEST01", "PAT02", "2012-06-19T00:00:00"
            ) %>% mutate(
              TRTEDTM = as_iso_dtm(TRTEDTM),
            )
            # nolint end

            actual_output <- derive_vars_trt_end(input,
                                                   input_ex,
                                                   new_var=TRTEDTM,
                                                   dtc=EXSTDTC,
                                                   date_imputation = "first",
                                                   time_imputation="first",
                                                   min_dates = NULL,
                                                   max_dates = NULL,
                                                   order = vars(EXSTDTC,EXSEQ),
                                                   subject_keys=vars(STUDYID,USUBJID),
                                                   filter_ex = NULL,
                                                   flag_imputation = "AUTO",
                                                   ord_filter="last")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

#test 3
test_that("Test3::Adds treatment start date - no imputations
          choosing the last ex date with filter of EXDOSE>0 applied,
          but new_var not named --DTM,
          expect warning thrown on name of new_var", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06-19", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDT,
              "TEST01", "PAT01", "2014-01-17T00:00:00",
              "TEST01", "PAT02", "2012-06-19T00:00:00"
            ) %>% mutate(
              TRTEDT = as_iso_dtm(TRTEDT),
            )
            # nolint end

            actual_output <- derive_vars_trt_end(input,
                                                   input_ex,
                                                   new_var=TRTEDT,
                                                   dtc=EXSTDTC,
                                                   date_imputation = "first",
                                                   time_imputation="first",
                                                   min_dates = NULL,
                                                   max_dates = NULL,
                                                   order = vars(EXSTDTC,EXSEQ),
                                                   subject_keys=vars(STUDYID,USUBJID),
                                                   #filter_ex = NULL,
                                                   flag_imputation = "AUTO",
                                                   filter_ex = EXDOSE>0 ,
                                                   ord_filter="last")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

#test 4
test_that("Test4::Adds treatment END date IMPUTING WHEN NECESSARY - and imputation flags
          choosing the last ex date with filter of EXDOSE>0 applied", {
            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01", 2 , 38,
              "TEST01", "PAT01", "2014-06-02", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDTM, ~TRTDTF,
              "TEST01", "PAT01", "2014-01-01T00:00:00", "D",
              "TEST01", "PAT02", "2012-06-19T00:00:00", NA_character_
            ) %>% mutate(
              TRTEDTM = as_iso_dtm(TRTEDTM),
            )
            # nolint end

            actual_output <- derive_vars_trt_end(input,
                                                   input_ex,
                                                   new_var=TRTEDTM,
                                                   dtc=EXSTDTC,
                                                   date_imputation = "first",
                                                   time_imputation="first",
                                                   min_dates = NULL,
                                                   max_dates = NULL,
                                                   order = vars(EXSTDTC,EXSEQ),
                                                   subject_keys=vars(STUDYID,USUBJID),
                                                   flag_imputation = "date",
                                                   filter_ex = EXDOSE>0 ,
                                                   ord_filter="LAST")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })
