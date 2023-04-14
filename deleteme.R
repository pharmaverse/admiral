library(admiral.test)
library(dplyr, warn.conflicts = FALSE)
n <- 3
VISITS = c("SCREENING", "WEEK 1", "WEEK 5", "WEEK 10")
vs_rows <- n * length(VISITS)

USUBJID <- rep(paste("ID", 1:n, sep =""), each = length(VISITS))
VISIT <- rep(VISITS, times = n)

# recreate vs with variable STUDYID, DOMAIN, USUBJID, VISIT, DOMAIN

vs <- tribble(
  ~USUBJID, ~VISIT,      ~TUSTRESC,
  "1",       "SCREENING", "TARGET",
  "1",       "WEEK 1",    "TARGET",
  "1",       "WEEK 5",    "TARGET",
  "2",       "SCREENING",    "NON-TARGET",
  "2",       "WEEK 1", "NON-TARGET",
  "2",       "WEEK 5", "NON-TARGET"
)

install.packages("datapasta")

admiral_vs
table(admiral_vs$STUDYID)
data("admiral_dm")

install.packages("datapasta")


# Merging all dm variables to vs
derive_vars_merged(
  admiral_vs,
  dataset_add = select(admiral_dm, -DOMAIN),
  by_vars = get_admiral_option("subject_keys")
) %>%
  select(STUDYID, USUBJID, VSTESTCD, VISIT, VSTPT, VSSTRESN, AGE, AGEU)
datapasta::dpasta(vs)
tibble::tribble(
  ~USUBJID,      ~VISIT,    ~TUSTRESC,
       "1", "SCREENING",     "TARGET",
       "1",    "WEEK 1",     "TARGET",
       "1",    "WEEK 5",     "TARGET",
       "2", "SCREENING", "NON-TARGET",
       "2",    "WEEK 1", "NON-TARGET",
       "2",    "WEEK 5", "NON-TARGET"
  )





# GET_ADMIRAL_OPTIONS EXAMPLE
n <- 2
data("admiral_vs")
data("admiral_dm")
set.seed(1235)

ids <- admiral_vs$USUBJID %>% unique() %>% sample(n)
visits <- {admiral_vs$VISIT %>% unique()}[c(3, 5)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[1:2]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]
dm <- admiral_dm %>% filter(USUBJID %in% ids) %>%
  select(STUDYID, DOMAIN, USUBJID, AGE, AGEU) %>%
  datapasta::dpasta()

tibble::tribble(
   ~STUDYID, ~DOMAIN,      ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01","DM", "01-701-1302",   61, "YEARS",
  "PILOT01",    "DM", "01-717-1344",   64, "YEARS"
  ) %>% datapasta::dpasta()


dm <- tibble::tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1302",   61, "YEARS",
  "PILOT01",    "DM", "17-1344",   64, "YEARS"
  )
vs <- admiral_vs %>%
  filter(USUBJID %in% ids,
         VISIT %in% visits,
         VSTESTCD %in% vstestcd,
         VSTPT %in% vstpt) %>%
  select(STUDYID, DOMAIN, USUBJID, VSTESTCD, VISIT, VSTPT, VSSTRESN) %>%
  datapasta::dpasta()

vs <- tibble::tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID, ~VSTESTCD,     ~VISIT,     ~VSTPT, ~VSSTRESN,
  "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE",    "LYING",        76,
  "PILOT01",    "VS", "01-1302",   "DIABP", "BASELINE", "STANDING",        87,
  "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2",    "LYING",        71,
  "PILOT01",    "VS", "01-1302",   "DIABP",   "WEEK 2", "STANDING",        79,
  "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE",    "LYING",        88,
  "PILOT01",    "VS", "17-1344",   "DIABP", "BASELINE", "STANDING",        86,
  "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2",    "LYING",        84,
  "PILOT01",    "VS", "17-1344",   "DIABP",   "WEEK 2", "STANDING",        82
  )

derive_vars_merged(
  vs,
  dataset_add = select(dm, -DOMAIN),
  by_vars = get_admiral_option("subject_keys")
)


# call derivationlibrary(dplyr, warn.conflicts = FALSE)
library(admiral.test)
data(admiral_ae)
data(admiral_adsl)
set.seed(235)
n <- 5

ids <- admiral_adsl$USUBJID %>% sample(n)
admiral_ae %>% select(
  STUDYID,
  DOMAIN,
  USUBJID,
  AESTDTC,
  AEENDTC
) %>%
  filter(
    USUBJID %in% ids
  )%>% datapasta::dpasta()

tibble::tribble(
                                        ~STUDYID, ~DOMAIN,      ~USUBJID,     ~AESTDTC,     ~AEENDTC,
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
                                  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
                                  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
                                  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
                                  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
                                  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
                                  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
                                  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
                                  ) %>% datapasta::dpasta()

adsl <- admiral_adsl %>% select(
  STUDYID,
  USUBJID,
  TRTSDT,
  TRTEDT
) %>%
  filter(
    USUBJID %in% ids
  ) %>%
  mutate(
    TRTSDT = as.character(TRTSDT),
    TRTEDT = as.character(TRTEDT)
  ) %>%
  datapasta::dpasta()


adsl <- tribble(
  ~STUDYID,  ~USUBJID,      ~TRTSDT,      ~TRTEDT,
 "PILOT01", "01-1307",           NA,           NA,
 "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
 "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
 "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
 "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
)  %>%
  mutate(
    across(TRTSDT:TRTEDT, as.Date)
  )

ae <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
)

adae <- ae %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT, TRTEDT),
    by_vars = exprs(USUBJID)
  )

## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
## one can add multiple variables in one go
call_derivation(
  dataset = adae,
  derivation = derive_vars_dt,
  variable_params = list(
    params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
    params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
  ),
  min_dates = exprs(TRTSDT),
  max_dates = exprs(TRTEDT)
)

## The above call using `call_derivation()` is equivalent to the following
adae %>%
  derive_vars_dt(
    new_vars_prefix = "AST",adsl <- tribble(
      ~STUDYID,  ~USUBJID,      ~TRTSDT,      ~TRTEDT,
      "PILOT01", "01-1307",           NA,           NA,
      "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
      "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
      "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
      "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
    )  %>%
      mutate(
        across(TRTSDT:TRTEDT, as.Date)
      ))

ae <- tribble(
      ~STUDYID, ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
      "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
      "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
      "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
      "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
      "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
      "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
      "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
    )

    adae <- ae %>%
      derive_vars_merged(
        dataset_add = adsl,
        new_vars = exprs(TRTSDT, TRTEDT),
        by_vars = exprs(USUBJID)
      )

    ## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
    ## one can add multiple variables in one go
    call_derivation(
      dataset = adae,
      derivation = derive_vars_dt,
      variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
      ),
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    )

    ## The above call using `call_derivation()` is equivalent to the following
    adae %>%
      derive_vars_dt(
        new_vars_prefix = "AST",
        dtc = AESTDTC,
        date_imputation = "first",
        min_dates = exprs(TRTSDT),
        max_dates = exprs(TRTEDT)
      ) %>%
      derive_vars_dt(
        new_vars_prefix = "AEN",
        dtc = AEENDTC,
        date_imputation = "last",
        min_dates = exprs(TRTSDT),
        max_dates = exprs(TRTEDT)
      )

    dtc = AESTDTC,
    date_imputation = "first",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  ) %>%
  derive_vars_dt(
    new_vars_prefix = "AEN",
    dtc = AEENDTC,
    date_imputation = "last",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  )






### params ---------------------------------------------------------------------

#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(admiral_adsl)
#'
adae <- admiral_ae[sample(1:nrow(admiral_ae), 1000), ] %>%
  select(USUBJID, AESTDTC, AEENDTC) %>%
  derive_vars_merged(
    dataset_add = admiral_adsl,
    new_vars = exprs(TRTSDT, TRTEDT),
    by_vars = exprs(USUBJID)
  )
#'
#' ## In order to derive both `ASTDT` and `AENDT` in `ADAE`, one can use `derive_vars_dt()`
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   )
#'
#' ## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
#' ## one can add multiple variables in one go.
#' ## The function arguments which are different from a variable to another (e.g. `new_vars_prefix`,
#' ## `dtc`, and `date_imputation`) are specified as a list of `params()` in the `variable_params`
#' ## argument of `call_derivation()`. All other arguments which are common to all variables
#' ## (e.g. `min_dates` and `max_dates`) are specified outside of `variable_params` (i.e. in `...`).
#' call_derivation(
#'   dataset = adae,
#'   derivation = derive_vars_dt,
#'   variable_params = list(
#'     params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
#'     params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
#'   ),
#'   min_dates = exprs(TRTSDT),
#'   max_dates = exprs(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the call using `derive_vars_dt()`
#' ## to derive variables `ASTDT` and `AENDT` separately at the beginning.
devtools::load_all()

library(dplyr, warn.conflicts = FALSE)

adsl <- tribble(
  ~STUDYID,  ~USUBJID,      ~TRTSDT,      ~TRTEDT,
  "PILOT01", "01-1307",           NA,           NA,
  "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
  "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
  "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
  "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
)  %>%
  mutate(
    across(TRTSDT:TRTEDT, as.Date)
  )



ae <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
  "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
)

adae <- ae %>%
  select(USUBJID, AESTDTC, AEENDTC) %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT, TRTEDT),
    by_vars = exprs(USUBJID)
  )

## In order to derive both `ASTDT` and `AENDT` in `ADAE`, one can use `derive_vars_dt()`
adae %>%
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    date_imputation = "first",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  ) %>%
  derive_vars_dt(
    new_vars_prefix = "AEN",
    dtc = AEENDTC,
    date_imputation = "last",
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  )


## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
## one can add multiple variables in one go.
## The function arguments which are different from a variable to another (e.g. `new_vars_prefix`,
## `dtc`, and `date_imputation`) are specified as a list of `params()` in the `variable_params`
## argument of `call_derivation()`. All other arguments which are common to all variables
## (e.g. `min_dates` and `max_dates`) are specified outside of `variable_params` (i.e. in `...`).
call_derivation(
  dataset = adae,
  derivation = derive_vars_dt,
  variable_params = list(
    params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
    params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
  ),
  min_dates = exprs(TRTSDT),
  max_dates = exprs(TRTEDT)
)

## The above call using `call_derivation()` is equivalent to the call using `derive_vars_dt()`
## to derive variables `ASTDT` and `AENDT` separately at the beginning.

# derive vars merged -----------------------------------------------------------
library(admiral.test)
library(dplyr, warn.conflicts = FALSE)
data("admiral_vs")
data("admiral_dm")
data("admiral_adsl")
data("admiral_ex")


n <- 2
set.seed(1235)
ids <- admiral_vs$USUBJID %>% unique() %>% sample(n)
visits <- {admiral_vs$VISIT %>% unique()}[c(3,  8, 14:16)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[c(2,6)]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]


ex <- admiral_ex %>%
  filter(
    USUBJID %in% ids
  ) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    EXSTDY,
    EXENDY,
    EXSTDTC,
    EXENDTC
  )%>% datapasta::dpasta()

adsl <- admiral_adsl %>%
  filter(
    USUBJID %in% ids
  ) %>%
  select(
    STUDYID,
    USUBJID,
    AGE,
    AGEU
  ) %>%  datapasta::dpasta()


vs <- admiral_vs %>% filter(
  USUBJID %in%  ids,
  #VISIT %in% visits,
  VSTESTCD %in% vstestcd#,
  #VSTPT %in% vstpt
)%>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VSTESTCD,
    VISIT,
    #    VSTPT,
    VSSTRESN,
    VSSTRESU,
    VSDTC
  ) %>% datapasta::dpasta()

dm <- admiral_dm %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    AGE,
    AGEU
  )%>% datapasta::dpasta()





## start here
library(dplyr, warn.conflicts = FALSE)
vs <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID, ~VSTESTCD,      ~VISIT, ~VSSTRESN, ~VSSTRESU,       ~VSDTC,
  "PILOT01",    "VS", "01-1302",  "HEIGHT", "SCREEN",     177.8,      "cm", "2013-08-20",
  "PILOT01",    "VS", "01-1302",  "WEIGHT", "SCREEN",     81.19,      "kg", "2013-08-20",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",  "BASE",      82.1,      "kg", "2013-08-29",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 2",     81.19,      "kg", "2013-09-15",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 4",     82.56,      "kg", "2013-09-24",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 6",     80.74,      "kg", "2013-10-08",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 8",      82.1,      "kg", "2013-10-22",
  "PILOT01",    "VS", "01-1302",  "WEIGHT",   "WEEK 12",      82.1,      "kg", "2013-11-05",
  "PILOT01",    "VS", "17-1344",  "HEIGHT", "SCREE",     163.5,      "cm", "2014-01-01",
  "PILOT01",    "VS", "17-1344",  "WEIGHT", "SCREE",     58.06,      "kg", "2014-01-01",
  "PILOT01",    "VS", "17-1344",  "WEIGHT",  "BASE",     58.06,      "kg", "2014-01-11",
  "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 2",     58.97,      "kg", "2014-01-24",
  "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 4",     57.97,      "kg", "2014-02-07",
  "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 6",     58.97,      "kg", "2014-02-19",
  "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 8",     57.79,      "kg", "2014-03-14"
)

dm <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1302",   61, "YEARS",
  "PILOT01",    "DM", "17-1344",   64, "YEARS"
)


# Merging all dm variables to vs
derive_vars_merged(
  vs,
  dataset_add = select(dm, -DOMAIN),
  by_vars = exprs(STUDYID, USUBJID)
) %>%
  select(STUDYID, USUBJID, VSTESTCD, VISIT, VSSTRESN, AGE, AGEU)


# Merge last weight to adsl
adsl <- tribble(
  ~STUDYID,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01", "01-1302",   61, "YEARS",
  "PILOT01", "17-1344",   64, "YEARS"
)


derive_vars_merged(
  adsl,
  dataset_add = vs,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(VSDTC),
  mode = "last",
  new_vars = exprs(LASTWGT = VSSTRESN, LASTWGTU = VSSTRESU),
  filter_add = VSTESTCD == "WEIGHT",
  match_flag = vsdatafl
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, LASTWGT, LASTWGTU, vsdatafl)


# Derive treatment start datetime (TRTSDTM)
ex <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID, ~EXSTDY, ~EXENDY,     ~EXSTDTC,     ~EXENDTC,
  "PILOT01",    "EX", "01-1302",       1,      18, "2013-08-29", "2013-09-15",
  "PILOT01",    "EX", "01-1302",      19,      69, "2013-09-16", "2013-11-05",
  "PILOT01",    "EX", "17-1344",       1,      14, "2014-01-11", "2014-01-24",
  "PILOT01",    "EX", "17-1344",      15,      63, "2014-01-25", "2014-03-14"
)
## Impute exposure start date to first date/time
ex_ext <- derive_vars_dtm(
  ex,
  dtc = EXSTDTC,
  new_vars_prefix = "EXST",
  highest_imputation = "M",
)
## Add first exposure datetime and imputation flags to adsl
derive_vars_merged(
  select(dm, STUDYID, USUBJID),
  dataset_add = ex_ext,
  by_vars = exprs(STUDYID, USUBJID),
  new_vars = exprs(TRTSDTM = EXSTDTM, TRTSDTF = EXSTDTF, TRTSTMF = EXSTTMF),
  order = exprs(EXSTDTM),
  mode = "first"
)


# Derive treatment end datetime (TRTEDTM)
## Impute exposure end datetime to last time, no date imputation
ex_ext <- derive_vars_dtm(
  ex,
  dtc = EXENDTC,
  new_vars_prefix = "EXEN",
  time_imputation = "last",
)
## Add last exposure datetime and imputation flag to adsl
derive_vars_merged(
  select(adsl, STUDYID, USUBJID),
  dataset_add = ex_ext,
  filter_add = !is.na(EXENDTM),
  by_vars = exprs(STUDYID, USUBJID),
  new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
  order = exprs(EXENDTM),
  mode = "last"
)

























# derive var merged cat---------------------------------------------------------

library(admiral.test)
library(dplyr, warn.conflicts = FALSE)
data("admiral_dm")
data("admiral_vs")



n <- 2
x <- sample(1:10000, size = 1)
set.seed(5584)

ids <- c(admiral_vs$USUBJID %>% unique() %>% sample(size = n), "01-701-1057")
visits <- {admiral_vs$VISIT %>% unique()}[c(1, 3:6, 14)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[c(2,6)]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]

vs <- admiral_vs %>% filter(
  USUBJID %in%  ids,
  VISIT %in% visits,
  VSTESTCD %in% vstestcd#,
  #VSTPT %in% vstpt
)%>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VISIT,
    VSTESTCD,
    VSSTRESN,

    VSSEQ,
    VSDTC
  ) %>% datapasta::dpasta()
dm <- admiral_dm %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    AGE,
    AGEU
  ) %>% datapasta::dpasta()




library(dplyr, warn.conflicts = FALSE)

vs <- tribble(
 ~STUDYID, ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSSEQ,       ~VSDTC,
"PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,     43, "2013-09-16",
"PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,    142, "2013-09-16",
"PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,    143, "2013-10-02",
"PILOT01",    "VS", "04-1127",    "WEEK 2",  "WEIGHT",     42.64,    144, "2013-10-16",
"PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,    145, "2013-10-30",
"PILOT01",    "VS", "04-1127",   "WEEK 26",  "WEIGHT",     43.09,    152, "2014-03-31",
"PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,     28, "2013-04-30",
"PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,     92, "2013-04-30",
"PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     93, "2013-05-14",
"PILOT01",    "VS", "06-1049",    "WEEK 2",  "WEIGHT",     58.29,     94, "2013-05-28",
"PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,     95, "2013-06-11"
)

dm <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1057",   59, "YEARS",
  "PILOT01",    "DM", "04-1127",   84, "YEARS",
  "PILOT01",    "DM", "06-1049",   60, "YEARS"
)
wgt_cat <- function(wgt) {
  case_when(
    wgt < 50 ~ "low",
    wgt > 90 ~ "high",
    TRUE ~ "normal"
  )
}

derive_var_merged_cat(
  dm,
  dataset_add = vs,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(VSDTC, VSSEQ),
  filter_add = VSTESTCD == "WEIGHT" & substr(VISIT, 1, 9) == "SCREENING",
  new_var = WGTBLCAT,
  source_var = VSSTRESN,
  cat_fun = wgt_cat,
  mode = "last"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, WGTBLCAT)



# defining a value for missing VS data
derive_var_merged_cat(
  dm,
  dataset_add = vs,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(VSDTC, VSSEQ),
  filter_add = VSTESTCD == "WEIGHT" & substr(VISIT, 1, 9) == "SCREENING",
  new_var = WGTBLCAT,
  source_var = VSSTRESN,
  cat_fun = wgt_cat,
  mode = "last",
  missing_value = "MISSING"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, WGTBLCAT)



# derive_var_merged_exist_flag--------------------------------------------------
library(admiral.test)
data("admiral_dm")
data("admiral_ae")
data("admiral_vs")
n <- 2
x <- sample(1:10000, size = 1)
set.seed(5584)

ids <- c(admiral_vs$USUBJID %>% unique() %>% sample(size = n), "01-701-1028")
visits <- {admiral_vs$VISIT %>% unique()}[c(1, 3, 6)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[c(2,6)]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]


dm <- admiral_dm %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    AGE,
    AGEU
  ) %>% datapasta::dpasta()
ae <- admiral_ae %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    AETERM,
    AEREL
  ) %>% datapasta::dpasta()

vs <- admiral_vs %>% filter(
  USUBJID %in%  ids,
  VISIT %in% visits,
  VSTESTCD %in% vstestcd#,
  #VSTPT %in% vstpt
)%>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VISIT,
    VSTESTCD,
    VSSTRESN,
    VSBLFL
  ) %>% datapasta::dpasta()

## start here
library(dplyr, warn.conflicts = FALSE)

dm <- tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM", "01-1028",   71, "YEARS",
  "PILOT01",    "DM", "04-1127",   84, "YEARS",
  "PILOT01",    "DM", "06-1049",   60, "YEARS"
)

ae <- tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID,    ~AETERM,     ~AEREL,
  "PILOT01",    "AE", "01-1028", "ERYTHEMA", "POSSIBLE",
  "PILOT01",    "AE", "01-1028", "PRURITUS", "PROBABLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "POSSIBLE",
  "PILOT01",    "AE", "06-1049",  "SYNCOPE", "PROBABLE"
)


derive_var_merged_exist_flag(
  dm,
  dataset_add = ae,
  by_vars = exprs(STUDYID, USUBJID),
  new_var = AERELFL,
  condition = AEREL == "PROBABLE"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, AERELFL)

vs <- tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSBLFL,
  "PILOT01",    "VS", "01-1028", "SCREENING",  "HEIGHT",     177.8,      NA,
  "PILOT01",    "VS", "01-1028", "SCREENING",  "WEIGHT",     98.88,      NA,
  "PILOT01",    "VS", "01-1028",  "BASELINE",  "WEIGHT",     99.34,     "Y",
  "PILOT01",    "VS", "01-1028",    "WEEK 4",  "WEIGHT",     98.88,      NA,
  "PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,      NA,
  "PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,      NA,
  "PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,     "Y",
  "PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,      NA,
  "PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,      NA,
  "PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,      NA,
  "PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     "Y",
  "PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,      NA
)
derive_var_merged_exist_flag(
  dm,
  dataset_add = vs,
  by_vars = exprs(STUDYID, USUBJID),
  filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
  new_var = WTBLHIFL,
  condition = VSSTRESN > 90,
  false_value = "N",
  missing_value = "M"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, WTBLHIFL)















# derive_var_merged_character --------------------------------------------------
ids <- c("01-701-1023",
        "01-701-1028",
        "01-701-1033")
data("admiral_dm")
data("admiral_ds")

dm <- admiral_dm %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    AGE,
    AGEU
  ) %>% datapasta::dpasta()

ds <- admiral_ds %>% filter(
  USUBJID %in% ids
) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VISIT,
    DSCAT,
    DSDECOD
  )%>% datapasta::dpasta()

library(dplyr, warn.conflicts = FALSE)
dm <- tribble(
   ~STUDYID, ~DOMAIN, ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01",    "DM",   "1023",   64, "YEARS",
  "PILOT01",    "DM",   "1028",   71, "YEARS",
  "PILOT01",    "DM",   "1033",   74, "YEARS"
)

ds <- tribble(
  ~STUDYID, ~DOMAIN, ~USUBJID,      ~VISIT,           ~DSCAT,          ~DSDECOD,
  "PILOT01",    "DS",   "1023",  "BASELINE", "PROT MILESTONE",      "RANDOMIZED",
  "PILOT01",    "DS",   "1023",    "WEEK 4",     "DISP EVENT",   "ADVERSE EVENT",
  "PILOT01",    "DS",   "1023",    "WEEK 4",    "OTHER EVENT", "FINAL LAB VISIT",
  "PILOT01",    "DS",   "1023", "RETRIEVAL",    "OTHER EVENT", "FINAL RET VISIT",
  "PILOT01",    "DS",   "1028",  "BASELINE", "PROT MILESTONE",      "RANDOMIZED",
  "PILOT01",    "DS",   "1028",   "WEEK 26",     "DISP EVENT",       "COMPLETED",
  "PILOT01",    "DS",   "1028",   "WEEK 26",    "OTHER EVENT", "FINAL LAB VISIT",
  "PILOT01",    "DS",   "1033",  "BASELINE", "PROT MILESTONE",      "RANDOMIZED",
  "PILOT01",    "DS",   "1033",    "WEEK 4",     "DISP EVENT", "TERM BY SPONSOR",
  "PILOT01",    "DS",   "1033",    "WEEK 4",    "OTHER EVENT", "FINAL LAB VISIT",
  "PILOT01",    "DS",   "1033", "RETRIEVAL",    "OTHER EVENT", "FINAL RET VISIT"
)

derive_var_merged_character(
  dm,
  dataset_add = ds,
  by_vars = exprs(STUDYID, USUBJID),
  new_var = DISPSTAT,
  filter_add = DSCAT == "DISP EVENT",
  source_var = DSDECOD,
  case = "title"
) %>%
  select(STUDYID, USUBJID, AGE, AGEU, DISPSTAT)













# derive_vars_merged_lookup ----------------------------------------------------
library(admiral.test)
library(dplyr, warn.conflicts = FALSE)
data("admiral_vs")
n <- 2
x <- sample(1:10000, size = 1)
set.seed(5584)

ids <- c(admiral_vs$USUBJID %>% unique() %>% sample(size = n), "01-701-1028")
visits <- {admiral_vs$VISIT %>% unique()}[c(1, 3, 6)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[c(2, 5, 6)]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]




vs <- admiral_vs %>%
  filter(
    USUBJID %in% ids,
    VISIT %in% visits,
    VSTESTCD %in% vstestcd#,
    #VSTPT %in% vstpt
  ) %>%
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VISIT,
    VSTESTCD,
    VSTEST
  ) %>% datapasta::dpasta()

library(dplyr, warn.conflicts = FALSE)
vs <- tribble(
   ~STUDYID, ~DOMAIN,  ~USUBJID,        ~VISIT, ~VSTESTCD,       ~VSTEST,
  "PILOT01",    "VS", "01-1028",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "01-1028",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "01-1028", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "01-1028",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "01-1028",      "WEEK 4",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "04-1325",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "04-1325", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "04-1325",      "WEEK 4",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",   "SCREENING",  "HEIGHT",      "Height",
  "PILOT01",    "VS", "10-1027",   "SCREENING",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027",    "BASELINE",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027",      "WEEK 4",    "TEMP", "Temperature",
  "PILOT01",    "VS", "10-1027", "SCREENING 1",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",    "BASELINE",  "WEIGHT",      "Weight",
  "PILOT01",    "VS", "10-1027",      "WEEK 4",  "WEIGHT",      "Weight"
)

param_lookup <- tribble(
 ~VSTESTCD,                   ~VSTEST, ~PARAMCD,                           ~PARAM,
   "SYSBP", "Systolic Blood Pressure",  "SYSBP", "Systolic Blood Pressure (mmHg)",
  "WEIGHT",                  "Weight", "WEIGHT",                    "Weight (kg)",
  "HEIGHT",                  "Height", "HEIGHT",                    "Height (cm)",
    "TEMP",             "Temperature",   "TEMP",                "Temperature (C)",
     "MAP",  "Mean Arterial Pressure",    "MAP",  "Mean Arterial Pressure (mmHg)",
     "BMI",         "Body Mass Index",    "BMI",        "Body Mass Index(kg/m^2)",
     "BSA",       "Body Surface Area",    "BSA",         "Body Surface Area(m^2)"
 )

derive_vars_merged_lookup(
  dataset = vs,
  dataset_add = param_lookup,
  by_vars = exprs(VSTESTCD),
  new_vars = exprs(PARAMCD, PARAM),
  print_not_mapped = TRUE
)




















# derive_var_anrind ------------------------------------------------------------
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(admiral.test)
data(admiral_vs)
n <- 2
x <- sample(1:10000, size = 1)
set.seed(5584)

ids <- c(admiral_vs$USUBJID %>% unique() %>% sample(size = n))
visits <- {admiral_vs$VISIT %>% unique()}[c(3, 6)]
vstestcd <- {admiral_vs$VSTESTCD %>% unique}[c(1, 3, 6)]
vstpt <- {admiral_vs$VSTPT %>% unique}[1:2]




vs <- admiral_vs %>%
  filter(
    USUBJID %in% ids,
    VISIT %in% visits,
    VSTESTCD %in% vstestcd#,
    #VSTPT %in% vstpt
  ) %>%
 select(
    STUDYID,
    DOMAIN,
    USUBJID,
    VISIT,
    VSTESTCD,
    VSSTRESN
  ) %>% datapasta::dpasta()

library(dplyr, warn.conflicts = FALSE)
vs <- tribble(
  ~STUDYID, ~DOMAIN,  ~USUBJID,     ~VISIT, ~VSTESTCD, ~VSSTRESN,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "DIABP",        82,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "DIABP",        82,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "DIABP",        86,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "DIABP",        80,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "DIABP",        88,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "DIABP",        84,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "PULSE",        68,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "PULSE",        66,
  "PILOT01",    "VS", "04-1025", "BASELINE",   "PULSE",        74,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "PULSE",        86,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "PULSE",        92,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",   "PULSE",        88,
  "PILOT01",    "VS", "04-1025", "BASELINE",  "WEIGHT",     55.52,
  "PILOT01",    "VS", "04-1025",   "WEEK 4",  "WEIGHT",     55.57,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "DIABP",        76,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "DIABP",        72,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "DIABP",        72,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "DIABP",        80,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "DIABP",        72,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "DIABP",        72,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "PULSE",        72,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "PULSE",        68,
  "PILOT01",    "VS", "11-1143", "BASELINE",   "PULSE",        68,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "PULSE",        60,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "PULSE",        80,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",   "PULSE",        80,
  "PILOT01",    "VS", "11-1143", "BASELINE",  "WEIGHT",     64.64,
  "PILOT01",    "VS", "11-1143",   "WEEK 4",  "WEIGHT",     66.45
) %>%
  mutate(
    PARAMCD = VSTESTCD,
    AVAL = VSSTRESN
  )

ref_ranges <- tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "DIABP",      60,     80,    40,    90,
  "PULSE",      60,    100,    40,   110
)

vs %>% filter(PARAMCD %in% c("PULSE", "DIABP")) %>%
  derive_vars_merged(ref_ranges, by_vars = exprs(PARAMCD)) %>%
  derive_var_anrind() %>%
  select(USUBJID, PARAMCD, AVAL, ANRLO:ANRIND)











































