# derive_vars_cat Test 1: Basic functionality with advs dataset

    Code
      result
    Output
      # A tibble: 10 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height  147. <160            2
       2 01-701-1023 Height  163. >=160           1
       3 01-701-1028 Height  178. >=160           1
       4 01-701-1033 Height  175. >=160           1
       5 01-701-1034 Height  155. <160            2
       6 01-701-1047 Height  149. <160            2
       7 01-701-1097 Height  169. >=160           1
       8 01-701-1111 Height  158. <160            2
       9 01-701-1115 Height  182. >=160           1
      10 01-701-1118 Height  180. >=160           1

# derive_vars_cat Test 2: Error when dataset is not a dataframe

    Argument `dataset` must be class <data.frame>, but is a list.

# derive_vars_cat Test 3: Error when definition is not an exprs object

    Must specify at least one column using the `~name` syntax.

# derive_vars_cat Test 4: Error when required columns are missing from dataset

    Required variable `VSTEST` is missing in `dataset`

# derive_vars_cat Test 5: Correct behavior when no conditions are met

    Code
      result
    Output
      # A tibble: 10 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height  147. <NA>           NA
       2 01-701-1023 Height  163. <NA>           NA
       3 01-701-1028 Height  178. <NA>           NA
       4 01-701-1033 Height  175. <NA>           NA
       5 01-701-1034 Height  155. <NA>           NA
       6 01-701-1047 Height  149. <NA>           NA
       7 01-701-1097 Height  169. <NA>           NA
       8 01-701-1111 Height  158. <NA>           NA
       9 01-701-1115 Height  182. <NA>           NA
      10 01-701-1118 Height  180. <NA>           NA

# derive_vars_cat Test 6: Overlapping conditions handled correctly

    Code
      result
    Output
      # A tibble: 10 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height  147. <NA>           NA
       2 01-701-1023 Height  163. >=160           1
       3 01-701-1028 Height  178. >=160           1
       4 01-701-1033 Height  175. >=160           1
       5 01-701-1034 Height  155. <NA>           NA
       6 01-701-1047 Height  149. <NA>           NA
       7 01-701-1097 Height  169. >=160           1
       8 01-701-1111 Height  158. >155            2
       9 01-701-1115 Height  182. >=160           1
      10 01-701-1118 Height  180. >=160           1

# derive_vars_cat Test 7: Handles missing values in dataset correctly

    Code
      result
    Output
      # A tibble: 10 x 5
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N
         <chr>       <chr>  <dbl> <chr>       <dbl>
       1 01-701-1015 Height   NA  <NA>           NA
       2 01-701-1023 Height   NA  <NA>           NA
       3 01-701-1028 Height   NA  <NA>           NA
       4 01-701-1033 Height   NA  <NA>           NA
       5 01-701-1034 Height   NA  <NA>           NA
       6 01-701-1047 Height  149. <160            2
       7 01-701-1097 Height  169. >=160           1
       8 01-701-1111 Height  158. <160            2
       9 01-701-1115 Height  182. >=160           1
      10 01-701-1118 Height  180. >=160           1

# derive_vars_cat Test 8: Error when condition is missing from `definition`

    Required variable `condition` is missing in `definition`

# derive_vars_cat Test 9: Conditions for multiple VSTESTs (Height and BILI)

    Code
      result
    Output
      # A tibble: 20 x 5
         USUBJID     VSTEST  AVAL AVALCAT1        AVALCA1N
         <chr>       <chr>  <dbl> <chr>              <dbl>
       1 01-701-1015 Height 147.  Height < 160           2
       2 01-701-1015 Weight  54.0 Weight < 66.68         4
       3 01-701-1023 Height 163.  Height >= 160          1
       4 01-701-1023 Weight  78.5 Weight >= 66.68        3
       5 01-701-1028 Height 178.  Height >= 160          1
       6 01-701-1028 Weight  98.9 Weight >= 66.68        3
       7 01-701-1033 Height 175.  Height >= 160          1
       8 01-701-1033 Weight  88.4 Weight >= 66.68        3
       9 01-701-1034 Height 155.  Height < 160           2
      10 01-701-1034 Weight  63.5 Weight < 66.68         4
      11 01-701-1047 Height 149.  Height < 160           2
      12 01-701-1047 Weight  66.2 Weight < 66.68         4
      13 01-701-1097 Height 169.  Height >= 160          1
      14 01-701-1097 Weight  78.0 Weight >= 66.68        3
      15 01-701-1111 Height 158.  Height < 160           2
      16 01-701-1111 Weight  60.3 Weight < 66.68         4
      17 01-701-1115 Height 182.  Height >= 160          1
      18 01-701-1115 Weight  78.7 Weight >= 66.68        3
      19 01-701-1118 Height 180.  Height >= 160          1
      20 01-701-1118 Weight  71.7 Weight >= 66.68        3

# derive_vars_cat Test 10: Adding an extra variable (flag) to the dataset

    Code
      result
    Output
      # A tibble: 10 x 6
         USUBJID     VSTEST  AVAL AVALCAT1 AVALCA1N extra_var
         <chr>       <chr>  <dbl> <chr>       <dbl> <lgl>    
       1 01-701-1015 Height  147. <160            2 FALSE    
       2 01-701-1023 Height  163. >=160           1 TRUE     
       3 01-701-1028 Height  178. >=160           1 TRUE     
       4 01-701-1033 Height  175. >=160           1 TRUE     
       5 01-701-1034 Height  155. <160            2 FALSE    
       6 01-701-1047 Height  149. <160            2 FALSE    
       7 01-701-1097 Height  169. >=160           1 TRUE     
       8 01-701-1111 Height  158. <160            2 FALSE    
       9 01-701-1115 Height  182. >=160           1 TRUE     
      10 01-701-1118 Height  180. >=160           1 TRUE     

