# consolidate_metadata Test 4: warn if variables differ

    Code
      consolidate_metadata(datasets = list(global = glob, study = stud), key_vars = exprs(
        id))
    Condition
      Warning:
      The variable names differ across the input datasets.
      i This message can be suppressed by setting `check_vars = "none"`.
    Output
      # A tibble: 3 x 4
        SOURCE    id val        var  
        <chr>  <dbl> <chr>      <chr>
      1 global     1 glob_val_1 <NA> 
      2 global     2 glob_val_2 <NA> 
      3 study      3 stud_val_3 abc  

