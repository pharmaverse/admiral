---
title: "Verify Templates Report"
format: html 
date: "2025-08-07"
#editor: visual
execute:
  keep-md: true
editor_options: 
  chunk_output_type: console
---


::: {.cell}

:::

::: {.cell}

:::



`` 11 `` Base ADaM datasets were found in the `admiral/inst/verify/old` directory. The datasets are \``` r input_dataset_names` ``



::: {.cell}
::: {.cell-output .cell-output-stdout}

```
## Verify Templates Check Complete! 
```


:::

::: {.cell-output .cell-output-stdout}

```
Date:  2025-08-07 
```


:::

::: {.cell-output .cell-output-stdout}

```
Run by:   
```


:::

::: {.cell-output .cell-output-stdout}

```
Git Ref:   
```


:::

::: {.cell-output .cell-output-stdout}

```
BASE:  Generated ADaM Datasets from Templates during Run 
```


:::

::: {.cell-output .cell-output-stdout}

```
COMPARE:  ADaM Datasets from pharmaverseadam  
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
There are columns in BASE that are not in COMPARE !!
```


:::

::: {.cell-output .cell-output-stdout}

````
<summary>❌ Dataset: adae</summary>

```

Differences found between the objects!

Summary of BASE and COMPARE
  ==================================================================
    PROPERTY             BASE                       COMP            
  ------------------------------------------------------------------
      Name         get(comp_dataset)          get(new_dataset)      
     Class     "tbl_df, tbl, data.frame"  "tbl_df, tbl, data.frame" 
    Rows(#)              1191                       1191            
   Columns(#)             107                        105            
  ------------------------------------------------------------------


There are columns in BASE that are not in COMPARE !!
  =========
   COLUMNS 
  ---------
   DOSEON  
    DOSEU  
  ---------

```

</details>

<summary>✅ Dataset: adcm</summary>

```

No issues were found!
```

</details>
````


:::

::: {.cell-output .cell-output-stderr}

```
Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
Not all Values Compared Equal
```


:::

::: {.cell-output .cell-output-stdout}

````
<summary>❌ Dataset: adeg</summary>

```

Differences found between the objects!

Summary of BASE and COMPARE
  ==================================================================
    PROPERTY             BASE                       COMP            
  ------------------------------------------------------------------
      Name         get(comp_dataset)          get(new_dataset)      
     Class     "tbl_df, tbl, data.frame"  "tbl_df, tbl, data.frame" 
    Rows(#)              78756                      78756           
   Columns(#)             108                        108            
  ------------------------------------------------------------------


Not all Values Compared Equal
  =============================
   Variable  No of Differences 
  -----------------------------
    EGSEQ          26717       
   EGSTRESU        24660       
    EGTPT          57540       
     ATPT          57540       
  -----------------------------


First 10 of 26717 rows are shown in table below
  ========================================
   VARIABLE  ..ROWNUMBER..  BASE  COMPARE 
  ----------------------------------------
    EGSEQ          1        127      1    
    EGSEQ          2        128      2    
    EGSEQ          3        129      3    
    EGSEQ          4        130      4    
    EGSEQ          5        131      5    
    EGSEQ          6        132      6    
    EGSEQ          7        133      7    
    EGSEQ          8        134      8    
    EGSEQ          9        135      9    
    EGSEQ         10        136     10    
  ----------------------------------------


First 10 of 24660 rows are shown in table below
  ===============================================
   VARIABLE  ..ROWNUMBER..    BASE      COMPARE  
  -----------------------------------------------
   EGSTRESU       12        BEATS/MIN  beats/min 
   EGSTRESU       13        BEATS/MIN  beats/min 
   EGSTRESU       14        BEATS/MIN  beats/min 
   EGSTRESU       16        BEATS/MIN  beats/min 
   EGSTRESU       17        BEATS/MIN  beats/min 
   EGSTRESU       18        BEATS/MIN  beats/min 
   EGSTRESU       20        BEATS/MIN  beats/min 
   EGSTRESU       21        BEATS/MIN  beats/min 
   EGSTRESU       22        BEATS/MIN  beats/min 
   EGSTRESU       24        BEATS/MIN  beats/min 
  -----------------------------------------------


First 10 of 57540 rows are shown in table below
  =================================================================
   VARIABLE  ..ROWNUMBER..  BASE              COMPARE              
  -----------------------------------------------------------------
    EGTPT         12         1    "AFTER LYING DOWN FOR 5 MINUTES" 
    EGTPT         13         2     "AFTER STANDING FOR 1 MINUTE"   
    EGTPT         14         3     "AFTER STANDING FOR 3 MINUTES"  
    EGTPT         16         1    "AFTER LYING DOWN FOR 5 MINUTES" 
    EGTPT         17         2     "AFTER STANDING FOR 1 MINUTE"   
    EGTPT         18         3     "AFTER STANDING FOR 3 MINUTES"  
    EGTPT         20         1    "AFTER LYING DOWN FOR 5 MINUTES" 
    EGTPT         21         2     "AFTER STANDING FOR 1 MINUTE"   
    EGTPT         22         3     "AFTER STANDING FOR 3 MINUTES"  
    EGTPT         24         1    "AFTER LYING DOWN FOR 5 MINUTES" 
  -----------------------------------------------------------------


First 10 of 57540 rows are shown in table below
  =================================================================
   VARIABLE  ..ROWNUMBER..  BASE              COMPARE              
  -----------------------------------------------------------------
     ATPT         12         1    "AFTER LYING DOWN FOR 5 MINUTES" 
     ATPT         13         2     "AFTER STANDING FOR 1 MINUTE"   
     ATPT         14         3     "AFTER STANDING FOR 3 MINUTES"  
     ATPT         16         1    "AFTER LYING DOWN FOR 5 MINUTES" 
     ATPT         17         2     "AFTER STANDING FOR 1 MINUTE"   
     ATPT         18         3     "AFTER STANDING FOR 3 MINUTES"  
     ATPT         20         1    "AFTER LYING DOWN FOR 5 MINUTES" 
     ATPT         21         2     "AFTER STANDING FOR 1 MINUTE"   
     ATPT         22         3     "AFTER STANDING FOR 3 MINUTES"  
     ATPT         24         1    "AFTER LYING DOWN FOR 5 MINUTES" 
  -----------------------------------------------------------------

```

</details>

<summary>✅ Dataset: adex</summary>

```

No issues were found!
```

</details>

<summary>✅ Dataset: adlb</summary>

```

No issues were found!
```

</details>

<summary>✅ Dataset: admh</summary>

```

No issues were found!
```

</details>

<summary>✅ Dataset: adpc</summary>

```

No issues were found!
```

</details>

<summary>✅ Dataset: adpp</summary>

```

No issues were found!
```

</details>

<summary>✅ Dataset: adppk</summary>

```

No issues were found!
```

</details>
````


:::

::: {.cell-output .cell-output-stderr}

```
Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
Not all Values Compared Equal
```


:::

::: {.cell-output .cell-output-stdout}

````
<summary>❌ Dataset: adsl</summary>

```

Differences found between the objects!

Summary of BASE and COMPARE
  ==================================================================
    PROPERTY             BASE                       COMP            
  ------------------------------------------------------------------
      Name         get(comp_dataset)          get(new_dataset)      
     Class     "tbl_df, tbl, data.frame"  "tbl_df, tbl, data.frame" 
    Rows(#)               306                        306            
   Columns(#)             54                         54             
  ------------------------------------------------------------------


Not all Values Compared Equal
  =============================
   Variable  No of Differences 
  -----------------------------
    SAFFL           52         
  -----------------------------


First 10 of 52 rows are shown in table below
  ========================================
   VARIABLE  ..ROWNUMBER..  BASE  COMPARE 
  ----------------------------------------
    SAFFL          7         N     <NA>   
    SAFFL         14         N     <NA>   
    SAFFL         18         N     <NA>   
    SAFFL         19         N     <NA>   
    SAFFL         28         N     <NA>   
    SAFFL         33         N     <NA>   
    SAFFL         38         N     <NA>   
    SAFFL         41         N     <NA>   
    SAFFL         43         N     <NA>   
    SAFFL         46         N     <NA>   
  ----------------------------------------

```

</details>

<summary>✅ Dataset: advs</summary>

```

No issues were found!
```

</details>
````


:::
:::
