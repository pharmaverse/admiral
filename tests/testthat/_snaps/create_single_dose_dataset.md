# create_single_dose_dataset Test 6: Error when a date variable contains NA values

    Code
      create_single_dose_dataset(input)
    Condition
      Error in `create_single_dose_dataset()`:
      ! The arguments `start_date` or `start_datetime` and `end_date` or `end_datetime` cannot contain NA values.
      i Please check "AENDT" for NA values.

# create_single_dose_dataset Test 7: Message for improper DT column names, ASTDT

    Code
      create_single_dose_dataset(input, start_date = ADTSTD, keep_source_vars = exprs(
        USUBJID, EXDOSFRQ, ADTSTD, ASTDTM, AENDT, AENDTM))
    Condition
      Error in `create_single_dose_dataset()`:
      ! The argument `start_date` is expected to have a name ending with "---DT".
      Please update as it does not follow the expected naming convention.

# create_single_dose_dataset Test 8: Message for improper DT column names, AENDT

    Code
      create_single_dose_dataset(input, end_date = ADTEND, )
    Condition
      Error in `create_single_dose_dataset()`:
      ! The argument `end_date` is expected to have a name like "xxxDT".
      Please check as it does not follow the expected naming convention.

# create_single_dose_dataset Test 9: error if no datetime and freq more than QD

    Code
      create_single_dose_dataset(input)
    Condition
      Error in `create_single_dose_dataset()`:
      ! There are dose frequencies more frequent than once a day, thus arguments `start_datetime` and `end_datetime` must be specified.

