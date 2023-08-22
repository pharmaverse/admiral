# signal_duplicate_records Test 2: dataset of duplicate records can be accessed using `get_duplicates_dataset()`

    Code
      get_duplicates_dataset()
    Output
      Duplicate records with respect to `USUBJID`.
      # A tibble: 4 x 3
        USUBJID COUNTRY  AAGE
      * <chr>   <chr>   <dbl>
      1 P01     GER        22
      2 P01     JPN        34
      3 P04     BRA        21
      4 P04     BRA        21

