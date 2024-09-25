# derive_vars_cat Test 1: Basic functionality with advs dataset

    Code
      result
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height 147.  <160            2
       2 01-701-1023 Height 163.  >=160           1
       3 01-701-1028 Height 178.  >=160           1
       4 01-701-1033 Height 175.  >=160           1
       5 01-701-1034 Height  NA   <NA>           NA
       6 01-701-1047 Height  NA   <NA>           NA
       7 01-701-1097 Height 169.  >=160           1
       8 01-701-1111 Height 158.  <160            2
       9 01-701-1115 Height 182.  >=160           1
      10 01-701-1118 Height 180.  >=160           1
      11 01-701-1015 Weight  54.0 <NA>           NA
      12 01-701-1023 Weight  78.5 <NA>           NA
      13 01-701-1028 Weight  98.9 <NA>           NA
      14 01-701-1033 Weight  88.4 <NA>           NA
      15 01-701-1034 Weight  NA   <NA>           NA
      16 01-701-1047 Weight  NA   <NA>           NA
      17 01-701-1097 Weight  78.0 <NA>           NA
      18 01-701-1111 Weight  60.3 <NA>           NA
      19 01-701-1115 Weight  78.7 <NA>           NA
      20 01-701-1118 Weight  71.7 <NA>           NA

---

    Code
      result2
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height 147.  <160            2
       2 01-701-1023 Height 163.  >=160           1
       3 01-701-1028 Height 178.  >=160           1
       4 01-701-1033 Height 175.  >=160           1
       5 01-701-1034 Height  NA   <NA>           NA
       6 01-701-1047 Height  NA   <NA>           NA
       7 01-701-1097 Height 169.  >=160           1
       8 01-701-1111 Height 158.  <160            2
       9 01-701-1115 Height 182.  >=160           1
      10 01-701-1118 Height 180.  >=160           1
      11 01-701-1015 Weight  54.0 <NA>           NA
      12 01-701-1023 Weight  78.5 <NA>           NA
      13 01-701-1028 Weight  98.9 <NA>           NA
      14 01-701-1033 Weight  88.4 <NA>           NA
      15 01-701-1034 Weight  NA   <NA>           NA
      16 01-701-1047 Weight  NA   <NA>           NA
      17 01-701-1097 Weight  78.0 <NA>           NA
      18 01-701-1111 Weight  60.3 <NA>           NA
      19 01-701-1115 Weight  78.7 <NA>           NA
      20 01-701-1118 Weight  71.7 <NA>           NA

# derive_vars_cat Test 2: Forgot to specify by_vars

    Column(s) in `definition` already exist in `dataset`.
    Did you forget to specify `by_vars`,
    or are you rerunning your code?

# derive_vars_cat Test 4: Error when definition is not an exprs object

    Argument `definition` must be a list of expressions but is a tibble.
    i To create a list of expressions use `exprs()`.

# derive_vars_cat Test 6: Correct behavior when no conditions are met

    Code
      result
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height 147.  <NA>           NA
       2 01-701-1023 Height 163.  <NA>           NA
       3 01-701-1028 Height 178.  <NA>           NA
       4 01-701-1033 Height 175.  <NA>           NA
       5 01-701-1034 Height  NA   <NA>           NA
       6 01-701-1047 Height  NA   <NA>           NA
       7 01-701-1097 Height 169.  <NA>           NA
       8 01-701-1111 Height 158.  <NA>           NA
       9 01-701-1115 Height 182.  <NA>           NA
      10 01-701-1118 Height 180.  <NA>           NA
      11 01-701-1015 Weight  54.0 <NA>           NA
      12 01-701-1023 Weight  78.5 <NA>           NA
      13 01-701-1028 Weight  98.9 <NA>           NA
      14 01-701-1033 Weight  88.4 <NA>           NA
      15 01-701-1034 Weight  NA   <NA>           NA
      16 01-701-1047 Weight  NA   <NA>           NA
      17 01-701-1097 Weight  78.0 <NA>           NA
      18 01-701-1111 Weight  60.3 <NA>           NA
      19 01-701-1115 Weight  78.7 <NA>           NA
      20 01-701-1118 Weight  71.7 <NA>           NA

# derive_vars_cat Test 7: Overlapping conditions handled correctly

    Code
      result
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height 147.  <160            3
       2 01-701-1023 Height 163.  <170            2
       3 01-701-1028 Height 178.  >=170           1
       4 01-701-1033 Height 175.  >=170           1
       5 01-701-1034 Height  NA   <NA>           NA
       6 01-701-1047 Height  NA   <NA>           NA
       7 01-701-1097 Height 169.  <170            2
       8 01-701-1111 Height 158.  <160            3
       9 01-701-1115 Height 182.  >=170           1
      10 01-701-1118 Height 180.  >=170           1
      11 01-701-1015 Weight  54.0 <NA>           NA
      12 01-701-1023 Weight  78.5 <NA>           NA
      13 01-701-1028 Weight  98.9 <NA>           NA
      14 01-701-1033 Weight  88.4 <NA>           NA
      15 01-701-1034 Weight  NA   <NA>           NA
      16 01-701-1047 Weight  NA   <NA>           NA
      17 01-701-1097 Weight  78.0 <NA>           NA
      18 01-701-1111 Weight  60.3 <NA>           NA
      19 01-701-1115 Weight  78.7 <NA>           NA
      20 01-701-1118 Weight  71.7 <NA>           NA

# derive_vars_cat Test 9: Conditions for multiple VSTESTs (Height and Weight)

    Code
      result
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1        AVALCA1N
         <chr>       <chr>  <dbl> <chr>              <dbl>
       1 01-701-1015 Height 147.  Height < 160           2
       2 01-701-1023 Height 163.  Height >= 160          1
       3 01-701-1028 Height 178.  Height >= 160          1
       4 01-701-1033 Height 175.  Height >= 160          1
       5 01-701-1034 Height  NA   <NA>                  NA
       6 01-701-1047 Height  NA   <NA>                  NA
       7 01-701-1097 Height 169.  Height >= 160          1
       8 01-701-1111 Height 158.  Height < 160           2
       9 01-701-1115 Height 182.  Height >= 160          1
      10 01-701-1118 Height 180.  Height >= 160          1
      11 01-701-1015 Weight  54.0 Weight < 66.68         2
      12 01-701-1023 Weight  78.5 Weight >= 66.68        1
      13 01-701-1028 Weight  98.9 Weight >= 66.68        1
      14 01-701-1033 Weight  88.4 Weight >= 66.68        1
      15 01-701-1034 Weight  NA   <NA>                  NA
      16 01-701-1047 Weight  NA   <NA>                  NA
      17 01-701-1097 Weight  78.0 Weight >= 66.68        1
      18 01-701-1111 Weight  60.3 Weight < 66.68         2
      19 01-701-1115 Weight  78.7 Weight >= 66.68        1
      20 01-701-1118 Weight  71.7 Weight >= 66.68        1

# derive_vars_cat Test 10: Adding an extra variable (flag) to the dataset

    Code
      result
    Output
      # A tibble: 20 x 6
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N extra_var
         <chr>       <chr>  <dbl> <chr>       <dbl> <lgl>    
       1 01-701-1015 Height 147.  <160            2 FALSE    
       2 01-701-1023 Height 163.  >=160           1 TRUE     
       3 01-701-1028 Height 178.  >=160           1 TRUE     
       4 01-701-1033 Height 175.  >=160           1 TRUE     
       5 01-701-1034 Height  NA   <NA>           NA NA       
       6 01-701-1047 Height  NA   <NA>           NA NA       
       7 01-701-1097 Height 169.  >=160           1 TRUE     
       8 01-701-1111 Height 158.  <160            2 FALSE    
       9 01-701-1115 Height 182.  >=160           1 TRUE     
      10 01-701-1118 Height 180.  >=160           1 TRUE     
      11 01-701-1015 Weight  54.0 <NA>           NA NA       
      12 01-701-1023 Weight  78.5 <NA>           NA NA       
      13 01-701-1028 Weight  98.9 <NA>           NA NA       
      14 01-701-1033 Weight  88.4 <NA>           NA NA       
      15 01-701-1034 Weight  NA   <NA>           NA NA       
      16 01-701-1047 Weight  NA   <NA>           NA NA       
      17 01-701-1097 Weight  78.0 <NA>           NA NA       
      18 01-701-1111 Weight  60.3 <NA>           NA NA       
      19 01-701-1115 Weight  78.7 <NA>           NA NA       
      20 01-701-1118 Weight  71.7 <NA>           NA NA       

# derive_vars_cat Test 12: definition has wrong shape

    Failed to convert `definition` to `tibble`. `definition` should be specified similarly to how you would specify a `tibble` using the `tribble()` function so it can be converted to `tibble` using `tribble()`.

