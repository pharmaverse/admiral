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











