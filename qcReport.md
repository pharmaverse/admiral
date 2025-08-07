# Verify Templates Report

2025-08-06

This is the QC report dated 2025-08-06.

The report compares the ADaM datasets in the `pharmaverseadam` package
with the ADaM `datasets generated from templates during this run`. The
datasets are compared using the `diffdf` package.

This run was initiated by on the Git ref.

    R version 4.5.1 (2025-06-13)
    Platform: x86_64-pc-linux-gnu
    Running under: Debian GNU/Linux 13 (trixie)

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.1 
    LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.1;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    time zone: America/Los_Angeles
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] stringr_1.5.1 tibble_3.3.0  diffdf_1.1.1 

    loaded via a namespace (and not attached):
     [1] digest_0.6.37     fastmap_1.2.0     xfun_0.52         magrittr_2.0.3   
     [5] glue_1.8.0        knitr_1.50        pkgconfig_2.0.3   htmltools_0.5.8.1
     [9] rmarkdown_2.29    lifecycle_1.0.4   cli_3.6.5         vctrs_0.6.5      
    [13] compiler_4.5.1    tools_4.5.1       evaluate_1.0.3    pillar_1.10.2    
    [17] yaml_2.3.10       rlang_1.1.6       jsonlite_2.0.0    stringi_1.8.7    

\<! – \#\| output: false \# why? –\>

`11` Base ADaM datasets were found in the `admiral/inst/verify/old`
directory. The datasets are \``` r input_dataset_names` ``

    ## Verify Templates Check Complete! 

    Date:  2025-08-06 

    Run by:   

    Git Ref:   

    BASE:  Generated ADaM Datasets from Templates during Run 

    COMPARE:  ADaM Datasets from pharmaverseadam  

    Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
    There are columns in BASE that are not in COMPARE !!

    <details>
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

    <details>
    <summary>✅ Dataset: adcm</summary>

    ```

    No issues were found!
    ```

    </details>

    Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
    Not all Values Compared Equal

    <details>
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

    <details>
    <summary>✅ Dataset: adex</summary>

    ```

    No issues were found!
    ```

    </details>

    <details>
    <summary>✅ Dataset: adlb</summary>

    ```

    No issues were found!
    ```

    </details>

    <details>
    <summary>✅ Dataset: admh</summary>

    ```

    No issues were found!
    ```

    </details>

    <details>
    <summary>✅ Dataset: adpc</summary>

    ```

    No issues were found!
    ```

    </details>

    <details>
    <summary>✅ Dataset: adpp</summary>

    ```

    No issues were found!
    ```

    </details>

    <details>
    <summary>✅ Dataset: adppk</summary>

    ```

    No issues were found!
    ```

    </details>

    Warning in diffdf(base = get(comp_dataset), compare = get(new_dataset)): 
    Not all Values Compared Equal

    <details>
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

    <details>
    <summary>✅ Dataset: advs</summary>

    ```

    No issues were found!
    ```

    </details>
