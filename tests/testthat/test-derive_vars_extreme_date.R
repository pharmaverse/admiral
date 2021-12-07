context("test-derive_vars_extreme_date")

#----test 1----
test_that("Test1::Adds treatment end date - time imputations
          choosing the LAST ex date with filter of EXDOSE>0 applied", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06-19", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDTM,~TRTETMF,
              "TEST01", "PAT01", "2014-01-17T23:59:59", "H",
              "TEST01", "PAT02", "2012-06-19T23:59:59",  "H"
            ) %>% mutate(
              TRTEDTM = as_iso_dtm(TRTEDTM),
            )
            # nolint end

            actual_output <- derive_vars_extreme_dtm(input,
                                                 input_ex,
                                                 new_var_prefix = "TRTE",
                                                 dtc=EXENDTC,
                                                 date_imputation = NULL,
                                                 time_imputation= "last",
                                                 min_dates = NULL,
                                                 max_dates = NULL,
                                                 order = vars(EXENDTC,EXSEQ),
                                                 subject_keys=vars(STUDYID,USUBJID),
                                                 flag_imputation = "AUTO",
                                                 filter_dataset_date = EXDOSE>0 ,
                                                 mode="LAST")

            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
           })

#----test 2----
 test_that("Test2::Adds treatment start date - time imputations
           choosing the LAST ex date with no filter applied
           -expect a warning thrown for no filter applied", {

            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01-02", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06-19", 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06-19", 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~TRTEDTM,~TRTETMF,
              "TEST01", "PAT01", "2014-06-19T23:59:59", "H",
              "TEST01", "PAT02", "2012-06-19T23:59:59",  "H"
            ) %>% mutate(
              TRTEDTM = as_iso_dtm(TRTEDTM),
            )
            # nolint end

            actual_output <- derive_vars_extreme_dtm(input,
                                                     input_ex,
                                                     new_var_prefix = "TRTE",
                                                     dtc=EXENDTC,
                                                     date_imputation = NULL,
                                                     time_imputation= "last",
                                                     min_dates = NULL,
                                                     max_dates = NULL,
                                                     order = vars(EXENDTC,EXSEQ),
                                                     subject_keys=vars(STUDYID,USUBJID),
                                                     flag_imputation = "AUTO",
                                                     filter_dataset_date = NULL ,
                                                     mode="LAST")

            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

#----test 3----
 test_that("Test3::Adds phase start date - with date & time imputations
          choosing the first ex date with filter of EXDOSE>0 applied", {
             input <- tibble::tribble(
               ~STUDYID, ~USUBJID,
               "TEST01", "PAT01",
              "TEST01", "PAT02"
             )

             input_ex <- tibble::tribble(
               ~STUDYID, ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE,
               "TEST01", "PAT01", "2014-01", 1 , 0,
               "TEST01", "PAT01", "2014-01-17", 2 , 38,
               "TEST01", "PAT01", "2014-06"   , 3 , 0,
               "TEST01", "PAT02", "2012-01-02", 1 , 54,
               "TEST01", "PAT02", "2012-01-17", 2 , 81,
               "TEST01", "PAT02", "2012-06"   , 3 , 54,
             )

             # nolint start
             expected_output <- tibble::tribble(
               ~STUDYID, ~USUBJID, ~APHO1STDTM , ~APHO1STDTF, ~APHO1STTMF,
               "TEST01", "PAT01", "2014-01-31T23:59:59", "D", "H",
               "TEST01", "PAT02", "2012-01-02T23:59:59", NA_character_, "H"
             ) %>% mutate(
               APHO1STDTM = as_iso_dtm(APHO1STDTM),
            )
             # nolint end

             actual_output <- derive_vars_extreme_dtm(input,
                                                      input_ex,
                                                      new_var_prefix = "APHO1ST",
                                                      dtc=EXENDTC,
                                                      date_imputation = "last",
                                                      time_imputation= "last",
                                                      min_dates = NULL,
                                                      max_dates = NULL,
                                                      order = vars(EXENDTC,EXSEQ),
                                                      subject_keys=vars(STUDYID,USUBJID),
                                                      flag_imputation = "AUTO",
                                                      filter_dataset_date = NULL ,
                                                      mode="first")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

###################REPEAT ALL THE ABOVE FOR DERIVE_VARS_EXTREME_DT

#----test 4----
 test_that("Test4::Adds treatment end date - NO imputations
          choosing the LAST ex date with filter of EXDOSE>0 applied", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE,
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
              "TEST01", "PAT01", "2014-01-17T23:59:59",
              "TEST01", "PAT02", "2012-06-19T23:59:59",
            ) %>% mutate(
              TRTEDT = as.Date(TRTEDT),
            )
            # nolint end

            actual_output <- derive_vars_extreme_dt(input,
                                                     input_ex,
                                                     new_var_prefix = "TRTE",
                                                     dtc=EXENDTC,
                                                     order = vars(EXENDTC,EXSEQ),
                                                     subject_keys=vars(STUDYID,USUBJID),
                                                     filter_dataset_date = EXDOSE>0 ,
                                                     mode="LAST")

            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

#----test 5----
 test_that("Test5::Adds phase start date - with date imputations
          choosing the first ex date with filter no filter applied", {
            input <- tibble::tribble(
              ~STUDYID, ~USUBJID,
              "TEST01", "PAT01",
              "TEST01", "PAT02"
            )

            input_ex <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE,
              "TEST01", "PAT01", "2014-01", 1 , 0,
              "TEST01", "PAT01", "2014-01-17", 2 , 38,
              "TEST01", "PAT01", "2014-06"   , 3 , 0,
              "TEST01", "PAT02", "2012-01-02", 1 , 54,
              "TEST01", "PAT02", "2012-01-17", 2 , 81,
              "TEST01", "PAT02", "2012-06"   , 3 , 54,
            )

            # nolint start
            expected_output <- tibble::tribble(
              ~STUDYID, ~USUBJID, ~APHO1STDT, ~APHO1STDTF,
              "TEST01", "PAT01", "2014-01-31T23:59:59", "D",
              "TEST01", "PAT02", "2012-01-02T23:59:59", NA_character_
            ) %>% mutate(
              APHO1STDT = as.Date(APHO1STDT)
            )
            # nolint end

            actual_output <- derive_vars_extreme_dt(input,
                                                     input_ex,
                                                     new_var_prefix = "APHO1ST",
                                                     dtc=EXENDTC,
                                                     date_imputation = "last",
                                                     min_dates = NULL,
                                                     max_dates = NULL,
                                                     order = vars(EXENDTC,EXSEQ),
                                                     subject_keys=vars(STUDYID,USUBJID),
                                                     flag_imputation = "AUTO",
                                                     filter_dataset_date = NULL ,
                                                     mode="first")


            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })
