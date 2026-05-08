# Load shared ATOX grading fixtures into the current test environment.
# This keeps large, repeated fixtures out of the test body while ensuring
# each `test_that()` block gets a fresh local copy.
local_atoxgr_test_data <- function(env = parent.frame()) {
  sys.source(
    testthat::test_path("testdata", "derive_var_atoxgr_data.R"),
    envir = env
  )
}

## Test 1: ATOXGR cannot be graded ----
test_that("derive_var_atoxgr Test 1: ATOXGR cannot be graded", {
  local_atoxgr_test_data()
  exp_out_1 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    NA_character_,      NA_character_,    NA_character_, NA_character_,  NA_character_,
    "Hypophosphatemia", NA_character_,    NA_character_, NA_character_,  NA_character_,
    NA_character_,      "Hyperglycemia",  NA_character_, NA_character_,  NA_character_,
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, NA_character_,  NA_character_,
    # ATOXGRL is ungradable so cannot say ATOXGR is normal
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",            NA_character_,
    # ATOXGRH is ungradable so cannot say ATOXGR is normal
    "Hypoglycemia",     "Hyperglycemia",  "0",           NA_character_,  NA_character_,
  )

  input_1 <- exp_out_1 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_1,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_1
  )
})

## Test 2: ATOXGR = 0 (normal) ----
test_that("derive_var_atoxgr Test 2: ATOXGR = 0 (normal)", {
  local_atoxgr_test_data()
  exp_out_2 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  "0",           "0",            "0",
    NA_character_,      "INR Increased",  NA_character_, "0",            "0",
    "Hypophosphatemia", NA_character_,    "0",           NA_character_,  "0",
  )

  input_2 <- exp_out_2 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_2,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_2
  )
})

## Test 3: ATOXGR > 0 (HYPER) ----
test_that("derive_var_atoxgr Test 3: ATOXGR > 0 (HYPER)", {
  local_atoxgr_test_data()
  exp_out_3 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, "1",            "1",
    "Hypoglycemia",     "Hyperglycemia",  "0",           "4",            "4",
    NA_character_,      "INR Increased",  NA_character_, "2",            "2",
  )

  input_3 <- exp_out_3 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_3,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_3
  )
})

## Test 4: ATOXGR < 0 (HYPO) ----
test_that("derive_var_atoxgr Test 4: ATOXGR < 0 (HYPO)", {
  local_atoxgr_test_data()
  exp_out_4 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  "3",           NA_character_,  "-3",
    "Hypoglycemia",     "Hyperglycemia",  "1",           "0",            "-1",
    "Hypophosphatemia", NA_character_,    "4",           NA_character_,  "-4",
  )

  input_4 <- exp_out_4 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_4,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_4
  )
})

# derive_var_atoxgr_dir

test_low <- function(expected, meta, high = "HIGH", low = "LOW") {
  input <- expected %>%
    select(-ATOXGRL)

  actual <- derive_var_atoxgr_dir(
    input,
    new_var = ATOXGRL,
    meta_criteria = meta,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    high_indicator = high,
    low_indicator = low,
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("TESTNUM")
  )
}

test_high <- function(expected, meta, high = "HIGH", low = "LOW") {
  input <- expected %>%
    select(-ATOXGRH)

  actual <- derive_var_atoxgr_dir(
    input,
    new_var = ATOXGRH,
    meta_criteria = meta,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = high,
    low_indicator = low,
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("TESTNUM")
  )
}

## Blood and lymphatic system disorders

### Anemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: <80 g/L
### Grade 2: <100 - 80g/L
### Grade 1: <LLN - 100 g/L


## Test 5a: CTCAEv4 Anemia ----
test_that("derive_var_atoxgr_dir Test 5a: CTCAEv4 Anemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_si, meta = atoxgr_criteria_ctcv4)
})

# CV_UNIT is g/dL so divide numeric values by 10

## Test 5b: CTCAEv4 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 5b: CTCAEv4 Anemia (CV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 5c: CTCAEv4 Anemia (SI + CV unit) ----
test_that("derive_var_atoxgr_dir Test 5c: CTCAEv4 Anemia (SI + CV unit)", {
  local_atoxgr_test_data()
  exp_anemia_sicv <- exp_anemia_cv %>%
    mutate(TESTNUM = TESTNUM + 16) %>%
    bind_rows(exp_anemia_si)

  atoxgr_criteria_ctcv4_siuscv <- atoxgr_criteria_ctcv4_uscv %>%
    bind_rows(atoxgr_criteria_ctcv4)

  test_low(expected = exp_anemia_sicv, meta = atoxgr_criteria_ctcv4_siuscv)
})

## Test 6a: CTCAEv5 Anemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 6a: CTCAEv5 Anemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_si, meta = atoxgr_criteria_ctcv5)
})

## Test 6b: CTCAEv5 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 6b: CTCAEv5 Anemia (CV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 6c: CTCAEv6 Anemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 6c: CTCAEv6 Anemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_si, meta = atoxgr_criteria_ctcv6)
})

## Test 6d: CTCAEv6 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 6d: CTCAEv6 Anemia (CV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_anemia_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Leukocytosis
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)

### Grade 3: >100,000/mm3


## Test 7a: CTCAEv4 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 7a: CTCAEv4 Leukocytosis (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_si, meta = atoxgr_criteria_ctcv4)
})


## Test 7b: CTCAEv4 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 7b: CTCAEv4 Leukocytosis (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 7c: CTCAEv4 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 7c: CTCAEv4 Leukocytosis (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv2, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 8a: CTCAEv5 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 8a: CTCAEv5 Leukocytosis (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_si, meta = atoxgr_criteria_ctcv5)
})

## Test 8b: CTCAEv5 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8b: CTCAEv5 Leukocytosis (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 8c: CTCAEv5 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8c: CTCAEv5 Leukocytosis (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv2, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 8d: CTCAEv6 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 8d: CTCAEv6 Leukocytosis (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_si, meta = atoxgr_criteria_ctcv6)
})

## Test 8e: CTCAEv6 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8e: CTCAEv6 Leukocytosis (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

## Test 8f: CTCAEv6 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8f: CTCAEv6 Leukocytosis (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_leuko_cv2, meta = atoxgr_criteria_ctcv6_uscv)
})


## Investigations

### Activated partial thromboplastin time prolonged
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: >2.5 x ULN
### Grade 2: >1.5 - 2.5 x ULN
### Grade 1: >ULN - 1.5 x ULN


## Test 9a: CTCAEv4 Act. partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 9a: CTCAEv4 Act. partial thromboplastin time prolonged", {
  local_atoxgr_test_data()
  test_high(expected = exp_aptt, meta = atoxgr_criteria_ctcv4)
})

## Test 9b: CTCAEv4 Activated partial thromboplastin time prolonged (SI + CV METDATA) ----
test_that(paste(
  "derive_var_atoxgr_dir Test 9b: CTCAEv4 Activated partial thromboplastin time prolonged",
  "(SI + CV METDATA)"
), {
  local_atoxgr_test_data()
  atoxgr_criteria_ctcv4_sicv <- atoxgr_criteria_ctcv4 %>%
    bind_rows(atoxgr_criteria_ctcv4_uscv)

  test_high(expected = exp_aptt, meta = atoxgr_criteria_ctcv4_sicv)
})

## Test 10a: CTCAEv5 Act. partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 10a: CTCAEv5 Act. partial thromboplastin time prolonged", {
  local_atoxgr_test_data()
  test_high(expected = exp_aptt, meta = atoxgr_criteria_ctcv5)
})

## Test 10b: CTCAEv6 Act. partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 10b: CTCAEv6 Act. partial thromboplastin time prolonged", {
  local_atoxgr_test_data()
  test_high(expected = exp_aptt, meta = atoxgr_criteria_ctcv6)
})

## Test 10c: CTCAEv6 Act. part. thromboplastin time prolonged (CV) ----
test_that("derive_var_atoxgr_dir Test 10c: CTCAEv6 Act. part. thromboplastin time prolonged (CV)", {
  local_atoxgr_test_data()
  test_high(expected = exp_aptt, meta = atoxgr_criteria_ctcv4_uscv)
})

### Alanine aminotransferase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN


## Test 11: CTCAEv4 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 11: CTCAEv4 Alanine aminotransferase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_alt_ctcv4, meta = atoxgr_criteria_ctcv4)
})

### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >3.0 - 5.0 x ULN if BL was normal OR >3.0 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 3.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal

## Test 12a: CTCAEv5 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 12a: CTCAEv5 Alanine aminotransferase increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_alt_ctcv5_norm <- exp_alt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_alt_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                            ~AVAL, ~BASE, ~ATOXGRH, ~TESTNUM,
    "Not a term",                            80,    40, NA,       13,
    NA_character_,                           60,    40, NA,       14,
    "Alanine aminotransferase Increased",   801,    40, "4",      15,
    "Alanine aminotransferase Increased",   800,    40, "3",      16,
    "Alanine aminotransferase Increased",   201,    40, "3",      17,
    "Alanine aminotransferase Increased",   200,    40, "2",      18,
    "Alanine aminotransferase Increased",   121,    40, "2",      19,
    "Alanine aminotransferase Increased",   120,    40, "1",      20,
    "Alanine aminotransferase Increased",    60,    40, "1",      21,
    "Alanine aminotransferase Increased",    59,    40, "0",      22,
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased",   100,    NA, NA,       23,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased",    NA,    40, NA,       24,
  ) %>%
    # set BASE to be HIGH (not normal) %>%
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_alt_ctcv5 <- exp_alt_ctcv5_norm %>%
    bind_rows(exp_alt_ctcv5_abn)

  test_high(expected = exp_alt_ctcv5, meta = atoxgr_criteria_ctcv5)
})

### NCICTCAEv6 same criteria as NCICTCAEv4 when BASELINE is normal

### Grade 4: >20.0 x ULN if baseline was normal or less than normal;
### >4.0 x baseline if baseline was >ULN
### Grade 3: >5.0 - 20.0 x ULN if baseline was normal or less than normal;
### >2.0 - 4.0 x baseline if baseline was >ULN up to 5 x ULN
### Grade 2: >3.0 - 5.0 x ULN if baseline was normal or less than normal;
### >1.5 - 2.0 x baseline if baseline was >ULN
### Grade 1: >ULN - 3.0 x ULN if BL was normal or less than normal;
### 1.0 - 1.5 x baseline if baseline was >ULN

## Test 12b: CTCAEv6 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 12b: CTCAEv6 Alanine aminotransferase increased", {
  local_atoxgr_test_data()
  # V6 and V5 criteria identical when BASELINE normal
  exp_alt_ctcv6_norm <- exp_alt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      AVALU = NA_character_,
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_alt_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,                            ~AVAL, ~BASE, ~ATOXGRH, ~TESTNUM,
    "Not a term",                            80,    40, NA,       13,
    NA_character_,                           60,    40, NA,       14,
    "Alanine aminotransferase Increased",   161,    40, "4",      15,
    "Alanine aminotransferase Increased",   160,    40, "3",      16,
    "Alanine aminotransferase Increased",    81,    40, "3",      17,
    "Alanine aminotransferase Increased",    80,    40, "2",      18,
    "Alanine aminotransferase Increased",    61,    40, "2",      19,
    "Alanine aminotransferase Increased",    60,    40, "1",      20,
    "Alanine aminotransferase Increased",    40,    40, "1",      21,
    "Alanine aminotransferase Increased",    39,    40, "0",      22,
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased",   100,    NA, NA,       23,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased",    NA,    40, NA,       24,
  ) %>%
    # set BASE to be HIGH (not normal)
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_alt_ctcv6 <- exp_alt_ctcv6_norm %>%
    bind_rows(exp_alt_ctcv6_abn)

  test_high(expected = exp_alt_ctcv6, meta = atoxgr_criteria_ctcv6)
})

### Alkaline phosphatase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN


## Test 13: CTCAEv4 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 13: CTCAEv4 Alkaline phosphatase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_alkp_ctcv4, meta = atoxgr_criteria_ctcv4)
})

### Alkaline phosphatase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >2.5 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR >2.0 - 2.5 x BL if BL was abnormal

## Test 14a: CTCAEv5 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 14a: CTCAEv5 Alkaline phosphatase increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_alkp_ctcv5_norm <- exp_alkp_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria

  exp_alkp_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                         ~AVAL,  ~BASE,  ~ATOXGRH, ~TESTNUM,
    "Not a term",                      80,     40,     NA,       13,
    NA_character_,                     60,     40,     NA,       14,
    "Alkaline phosphatase increased",  801,    40,     "4",      15,
    "Alkaline phosphatase increased",  800,    40,     "3",      16,
    "Alkaline phosphatase increased",  201,    40,     "3",      17,
    "Alkaline phosphatase increased",  200,    40,     "2",      18,
    "Alkaline phosphatase increased",  101,    40,     "2",      19,
    "Alkaline phosphatase increased",  100,    40,     "1",      20,
    "Alkaline phosphatase increased",  80,     40,     "1",      21,
    "Alkaline phosphatase increased",  79,     40,     "0",      22,
    # ANRHI missing - cannot grade
    "Alkaline phosphatase increased",  100,    NA,     NA,       23,
    # AVAL missing cannot grade
    "Alkaline phosphatase increased",  NA,     40,     NA,       24,
  ) %>%
    # set BASE to be abnormal (HIGH) and create FLAG
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_alkp_ctcv5 <- exp_alkp_ctcv5_norm %>%
    bind_rows(exp_alkp_ctcv5_abn)

  test_high(expected = exp_alkp_ctcv5, meta = atoxgr_criteria_ctcv5)
})

### Alkaline phosphatase increased
### NCICTCAEv6 different to NCICTCAEv5
### Grade 1: >Baseline and ULN


## Test 14b: CTCAEv6 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 14b: CTCAEv6 Alkaline phosphatase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_alkp_ctcv6, meta = atoxgr_criteria_ctcv6)
})

## Test 14c: CTCAEv6 Alkaline phosphatase increased (CV unit) ----
test_that("derive_var_atoxgr_dir Test 14c: CTCAEv6 Alkaline phosphatase increased (CV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_alkp_ctcv6, meta = atoxgr_criteria_ctcv6_uscv)
})


### Aspartate aminotransferase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN


## Test 15: CTCAEv4 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 15: CTCAEv4 Aspartate aminotransferase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_ast_ctcv4, meta = atoxgr_criteria_ctcv4)
})

### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >3.0 - 5.0 x ULN if BL was normal OR >3.0 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 3.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal

## Test 16a: CTCAEv5 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 16a: CTCAEv5 Aspartate aminotransferase increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_ast_ctcv5_norm <- exp_ast_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_ast_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                              ~AVAL,  ~BASE,  ~ATOXGRH, ~TESTNUM,
    "Not a term",                           80,     40,     NA,       13,
    NA_character_,                          60,     40,     NA,       14,
    "Aspartate aminotransferase Increased", 801,    40,     "4",      15,
    "Aspartate aminotransferase Increased", 800,    40,     "3",      16,
    "Aspartate aminotransferase Increased", 201,    40,     "3",      17,
    "Aspartate aminotransferase Increased", 200,    40,     "2",      18,
    "Aspartate aminotransferase Increased", 121,    40,     "2",      19,
    "Aspartate aminotransferase Increased", 120,    40,     "1",      20,
    "Aspartate aminotransferase Increased", 60,     40,     "1",      21,
    "Aspartate aminotransferase Increased", 59,     40,     "0",      22,
    # ANRHI missing - cannot grade
    "Aspartate aminotransferase Increased", 100,    NA,     NA,       23,
    # AVAL missing cannot grade
    "Aspartate aminotransferase Increased", NA,     40,     NA,       24,
  ) %>%
    # set BASE to be abnormal (HIGH or H)  and create FLAG
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_ast_ctcv5 <- exp_ast_ctcv5_norm %>%
    bind_rows(exp_ast_ctcv5_abn)

  test_high(expected = exp_ast_ctcv5, meta = atoxgr_criteria_ctcv5)
})

### NCICTCAEv6 same criteria as NCICTCAEv5 when BASELINE is normal
### Grade 4: >20.0 x ULN if baseline was normal or less than normal;
### >4.0 x baseline if baseline was >ULN"
### Grade 3: >5.0 - 20.0 x ULN if baseline was normal or less than normal;
### >2.0 - 4.0 x baseline if baseline was >ULN up to 5 x ULN
### Grade 2: >3.0 - 5.0 x ULN if baseline was normal or less than normal;
### >1.5 - 2.0 x baseline if baseline was >ULN
### Grade 1: >ULN - 3.0 x ULN if baseline was normal or less than normal;
### 1.0 - 1.5 x baseline if baseline was >ULN

## Test 16b: CTCAEv6 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 16b: CTCAEv6 Aspartate aminotransferase increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_ast_ctcv6_norm <- exp_ast_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_ast_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,                              ~AVAL,  ~BASE,  ~ATOXGRH, ~TESTNUM,
    "Not a term",                           80,     40,     NA,       13,
    NA_character_,                          60,     40,     NA,       14,
    "Aspartate aminotransferase Increased", 161,    40,     "4",      15,
    "Aspartate aminotransferase Increased", 160,    40,     "3",      16,
    "Aspartate aminotransferase Increased", 81,     40,     "3",      17,
    "Aspartate aminotransferase Increased", 80,     40,     "2",      18,
    "Aspartate aminotransferase Increased", 61,     40,     "2",      19,
    "Aspartate aminotransferase Increased", 60,     40,     "1",      20,
    "Aspartate aminotransferase Increased", 40,     40,     "1",      21,
    "Aspartate aminotransferase Increased", 39,     40,     "0",      22,
    # ANRHI missing - cannot grade
    "Aspartate aminotransferase Increased", 100,    NA,     NA,       23,
    # AVAL missing cannot grade
    "Aspartate aminotransferase Increased", NA,     40,     NA,       24,
  ) %>%
    # set BASE to be abnormal (HIGH or H)  and create FLAG
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_ast_ctcv6 <- exp_ast_ctcv6_norm %>%
    bind_rows(exp_ast_ctcv6_abn)

  test_high(expected = exp_ast_ctcv6, meta = atoxgr_criteria_ctcv6)
})

### Blood bilirubin increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >10.0 x ULN
### Grade 3: >3.0 - 10.0 x ULN
### Grade 2: >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN


## Test 17: CTCAEv4 Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 17: CTCAEv4 Blood bilirubin increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_bili_ctcv4, meta = atoxgr_criteria_ctcv4)
})

## Test 18a: CTCAEv5  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 18a: CTCAEv5  Blood bilirubin increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_bili_ctcv5_norm <- exp_bili_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL",
      TESTNUM = TESTNUM + 12
    )

  # create records with abnormal BASE then add records with normal BASE
  exp_bili_ctcv5 <- exp_bili_ctcv4 %>%
    # set BASE to ANRHI then make ANRHI < BASE
    mutate(
      BASE = ANRHI,
      ANRHI = ANRHI - 1,
      BNRIND = "HIGH"
    ) %>%
    bind_rows(exp_bili_ctcv5_norm)

  test_high(expected = exp_bili_ctcv5, meta = atoxgr_criteria_ctcv5)
})

### NCICTCAEv6 small difference to NCICTCAEv5 when BASELINE is abnormal
### Blood bilirubin increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >10.0 x ULN if baseline was normal or less than normal;
### >10.0 x baseline if baseline was >ULN
### Grade 3: >3.0 - 10.0 x ULN if baseline was normal or less than normal;
### >2.5 - 10.0 x baseline if baseline was >ULN
### Grade 2: >1.5 - 3.0 x ULN if baseline was normal or less than normal;
### >1.5 - 2.5 x baseline if baseline was >ULN
### Grade 1: >ULN - 1.5 x ULN if baseline was normal or less than normal;
### 1.0 - 1.5 x baseline if baseline was >ULN

# V5 and V6 criteria identical when BASELINE normal

# create records with abnormal BASE then add records with normal BASE

## Test 18b: CTCAEv6  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 18b: CTCAEv6  Blood bilirubin increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_bili_ctcv6, meta = atoxgr_criteria_ctcv6)
})

### CD4 Lymphocytes decreased
### NCICTCAEv5 and NCICTCAEv6 same criteria as NCICTCAEv4
### SI unit (10^9/L)
### Grade 4: <0.05 x 10e9 /L
### Grade 3: <0.2 - 0.05 x 10e9 /L
### Grade 2: <0.5 - 0.2 x 10e9 /L
### Grade 1: <LLN - 0.5 x 10e9 /L
### USCV unit (1/uL = /mm3)
### Grade 4: <50/mm3
### Grade 3: <200 x 50/mm3
### Grade 2: <500 - 200/mm3
### Grade 1: <LLN - 500/mm3



## Test 19a: CTCAEv4 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 19a: CTCAEv4 CD4 Lymphocytes decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_si, meta = atoxgr_criteria_ctcv4)
})

## Test 19b: CTCAEv4 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 19b: CTCAEv4 CD4 Lymphocytes decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 20a: CTCAEv5 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 20a: CTCAEv5 CD4 Lymphocytes decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_si, meta = atoxgr_criteria_ctcv5)
})

## Test 20b: CTCAEv5 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 20b: CTCAEv5 CD4 Lymphocytes decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 20c: CTCAEv6 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 20c: CTCAEv6 CD4 Lymphocytes decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_si, meta = atoxgr_criteria_ctcv6)
})

## Test 20d: CTCAEv6 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 20d: CTCAEv6 CD4 Lymphocytes decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Cholesterol high
### NCICTCAEv5 same criteria as NCICTCAEv4
### SI unit (mmol/L)
### Grade 4: >12.92 mmol/L
### Grade 3: >10.34 - 12.92 mmol/L
### Grade 2: >7.75 - 10.34 mmol/L
### Grade 1: >ULN - 7.75 mmol/L


### SI unit (mg/dL)
### Grade 4: >500 mg/dL
### Grade 3: >400 - 500 mg/dL
### Grade 2: >300 - 400 mg/dL
### Grade 1: >ULN - 300 mg/dL



## Test 21a: CTCAEv4 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 21a: CTCAEv4 Cholesterol high (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_si, meta = atoxgr_criteria_ctcv4)
})

## Test 21b: CTCAEv4 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 21b: CTCAEv4 Cholesterol high (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 22a: CTCAEv5 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 22a: CTCAEv5 Cholesterol high (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_si, meta = atoxgr_criteria_ctcv5)
})

## Test 22b: CTCAEv5 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 22b: CTCAEv5 Cholesterol high (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 22c: CTCAEv6 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 22c: CTCAEv6 Cholesterol high (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_si, meta = atoxgr_criteria_ctcv6)
})

## Test 22d: CTCAEv6 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 22d: CTCAEv6 Cholesterol high (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_choles_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### CPK increased
### NCICTCAEv5 same criteria as NCICTCAEv4
### Grade 4: >10.0 x ULN
### Grade 3: >5.0 - 10.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN


## Test 23: CTCAEv4 CPK increased ----
test_that("derive_var_atoxgr_dir Test 23: CTCAEv4 CPK increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_cpk, meta = atoxgr_criteria_ctcv4)
})

## Test 24: CTCAEv5 CPK increased ----
test_that("derive_var_atoxgr_dir Test 24: CTCAEv5 CPK increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_cpk, meta = atoxgr_criteria_ctcv5)
})

### Creatinine increased (NCICTCv4)
### NCICTCAEv5 same criteria as NCICTCAEv4 except for Grade 1
### Grade 4: >6.0 x ULN
### Grade 3: >3.0 baseline; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x baseline; >1.5 - 3.0 x ULN
### Grade 1: >1 - 1.5 x baseline; >ULN - 1.5 x ULN

# create flag to remove obs not relevant for NCI-CTCAEv5

## Test 25: CTCAEv4 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 25: CTCAEv4 Creatinine increased", {
  local_atoxgr_test_data()
  exp_creatn <- exp_creatn %>%
    filter(V4 == "Y")

  test_high(expected = exp_creatn, meta = atoxgr_criteria_ctcv4)
})

### Creatinine increased (NCICTCv5)
### NCICTCAEv5 same criteria as NCICTCAEv4 except for Grade 1
### Grade 4: >6.0 x ULN
### Grade 3: >3.0 baseline; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x baseline; >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

## Test 26a: CTCAEv5 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 26a: CTCAEv5 Creatinine increased", {
  local_atoxgr_test_data()
  exp_creatn <- exp_creatn %>%
    filter(V5 == "Y")

  test_high(expected = exp_creatn, meta = atoxgr_criteria_ctcv5)
})

### Creatinine increased (NCICTCv6)

### Grade 4: >6.0 x ULN
### Grade 3: >3.0 x if baseline is below LLN; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x if baseline is below LLN; >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

## Test 26b: CTCAEv6 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 26b: CTCAEv6 Creatinine increased", {
  local_atoxgr_test_data()
  exp_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~BNRIND,        ~ATOXGRH,      ~TESTNUM,
    "Not a term",            80,     80,    40,     "NORMAL",       NA_character_, 1,
    NA_character_,           60,     60,    40,     "NORMAL",       NA_character_, 2,
    # GRADE derived from AVAL against ANRHI
    "Creatinine increased",  241,    40,    40,     "NORMAL",       "4",           3,
    "Creatinine increased",  240,    40,    40,     "NORMAL",       "3",           4,
    "Creatinine increased",  121,    40,    40,     "NORMAL",       "3",           5,
    "Creatinine increased",  120,    40,    40,     "NORMAL",       "2",           6,
    "Creatinine increased",  61,     40,    40,     "NORMAL",       "2",           7,
    "Creatinine increased",  60,     40,    40,     "NORMAL",       "1",           8,
    "Creatinine increased",  41,     40,    40,     "NORMAL",       "1",           9,
    "Creatinine increased",  40,     40,    40,     "NORMAL",       "0",           10,
    # GRADE derived from AVAL against BASE when baseline is below LLN
    "Creatinine increased",  121,    40,    41,     "LOW",          "3",           11,
    "Creatinine increased",  120,    40,    41,     "LOW",          "2",           12,
    "Creatinine increased",  61,     40,    41,     "LOW",          "2",           13,
    "Creatinine increased",  60,     40,    60,     "LOW",          "0",           14,
    # BASE missing - AVAL <= ANRLO can grade as NORMAL
    "Creatinine increased",  42,     NA,    42,     NA_character_,  "0",           15,
    # ANRHI missing - baseline is below LLN and AVAL <= 1.5 x BASE but cannot grade as NORMAL
    "Creatinine increased",  60,     40,    NA,     "LOW",          NA_character_, 16,
    # AVAL missing cannot grade
    "Creatinine increased",  NA,     0,     40,     "NORMAL",       NA_character_, 17,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_creatn, meta = atoxgr_criteria_ctcv6)
})

### Fibrinogen decreased (NCICTCv4)
### Grade 4: <0.25 x LLN or 75% decrease from baseline or absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN or 50 - <75% decrease from baseline
### Grade 2: <0.75 - 0.5 x LLN or 25 - <50% decrease from baseline
### Grade 1: <1.0 - 0.75 x LLN or <25% decrease from baseline

## Test 27: CTCAEv4 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 27: CTCAEv4 Fibrinogen decreased", {
  local_atoxgr_test_data()
  exp_fib <- tibble::tribble(
    ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~PCHG,  ~AVALU,  ~ATOXGRL, ~TESTNUM,
    "Not a term",            9,      10,     40,     "g/L",   NA,       1,
    NA_character_,           10,     10,     40,     "g/L",   NA,       2,
    # Satisfies < 0.5 for grade 4 - other criteria missing
    "Fibrinogen decreased",  0.49,   NA,     NA,     "g/L",   "4",      3,
    # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
    "Fibrinogen decreased",  0.49,   1,      -51,    "g/L",   "4",      4,
    # Satisfies < 0.25*LLN for grade 4 - PCHG missing
    "Fibrinogen decreased",  0.5,    2.1,    NA,     "g/L",   "4",      5,
    # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
    "Fibrinogen decreased",  0.5,    2.1,    -51,    "g/L",   "4",      6,
    # Satisfies <=75% decrease for grade 4 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -75,    "g/L",   "4",      7,
    # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
    "Fibrinogen decreased",  1,      0.49,   -75,    "g/L",   "4",      8,
    # Satisfies < 0.5*LLN for grade 3 - PCHG missing
    "Fibrinogen decreased",  1,      2.1,    NA,     "g/L",   "3",      9,
    # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
    "Fibrinogen decreased",  1,      2.1,    -49,    "g/L",   "3",      10,
    # Satisfies <=50% decrease for grade 3 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -50,    "g/L",   "3",      11,
    # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
    "Fibrinogen decreased",  1,      2,      -50,    "g/L",   "3",      12,
    # Satisfies < 0.75*LLN for grade 2 - PCHG missing
    "Fibrinogen decreased",  1.5,    2.1,    NA,     "g/L",   "2",      13,
    # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
    "Fibrinogen decreased",  1.5,    2.1,    -10,    "g/L",   "2",      14,
    # Satisfies <=25% for grade 2 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -25,    "g/L",   "2",      15,
    # Satisfies <=25% for grade 2 - LLN satisfies grade 1
    "Fibrinogen decreased",  1.5,    1.6,    -25,    "g/L",   "2",      16,
    # Satisfies < LLN for grade 1 - PCHG missing
    "Fibrinogen decreased",  2,      2.1,    NA,     "g/L",   "1",      17,
    # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
    "Fibrinogen decreased",  2,      2.1,    10,     "g/L",   "1",      18,
    # Satisfies % decrease for grade 1 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -1,     "g/L",   "1",      19,
    # Satisfies % decrease for grade 1 - AVAL = LLN
    "Fibrinogen decreased",  1.5,    1.5,    -1,     "g/L",   "1",      20,
    # Satisfies grade 0 - AVAL >= LLN AND no % descrease
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/L",   "0",      21,
    # AVAL >= LLN BUT PCT missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    1.5,    NA,     "g/L",   NA,       22,
    # PCT >= 0 BUT LLN missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    NA,     10,     "g/L",   NA,       23,
    # AVAL missing cannot grade
    "Fibrinogen decreased",  NA,     1.5,    10,     "g/L",   NA,       24,
    # wrong unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/dL",  NA,       25,
    # missing unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      NA,      NA,       26,
  )

  test_low(expected = exp_fib, meta = atoxgr_criteria_ctcv4)
})

### Fibrinogen decreased (NCICTCv5 + NCICTCv6)
### Grade 4: <0.25 x LLN OR if abnormal, 75% dec. from BL OR absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN OR if abnormal, 50 - <75% dec. from BL
### Grade 2: <0.75 - 0.5 x LLN OR if abnormal, 25 - <50% dec. from BL
### Grade 1: <1.0 - 0.75 x LLN OR if abnormal, <25% dec. from BL


## Test 28a: CTCAEv5 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 28a: CTCAEv5 Fibrinogen decreased", {
  local_atoxgr_test_data()
  test_low(expected = exp_fib, meta = atoxgr_criteria_ctcv5)
})

## Test 28b: CTCAEv6 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 28b: CTCAEv6 Fibrinogen decreased", {
  local_atoxgr_test_data()
  test_low(expected = exp_fib, meta = atoxgr_criteria_ctcv6)
})

### GGT increased (NCICTCv4)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN


## Test 29: CTCAEv4 GGT increased ----
test_that("derive_var_atoxgr_dir Test 29: CTCAEv4 GGT increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_ggt_ctcv4, meta = atoxgr_criteria_ctcv4)
})

### GGT increased (NCICTCv5)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >2.5 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR >2.0 - 2.5 x BL if BL was abnormal

## Test 30a: CTCAEv5 GGT increased ----
test_that("derive_var_atoxgr_dir Test 30a: CTCAEv5 GGT increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_ggt_ctcv5_norm <- exp_ggt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_ggt_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRLO, ~BASE, ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    0,      40,    NA,       13,
    NA_character_,   60,    0,      40,    NA,       14,
    "GGT increased", 801,   0,      40,    "4",      15,
    "GGT increased", 800,   0,      40,    "3",      16,
    "GGT increased", 201,   0,      40,    "3",      17,
    "GGT increased", 200,   0,      40,    "2",      18,
    "GGT increased", 101,   0,      40,    "2",      19,
    "GGT increased", 100,   0,      40,    "1",      20,
    "GGT increased", 81,    0,      40,    "1",      21,
    "GGT increased", 80,    0,      40,    "0",      22,
    # ANRHI missing - cannot grade
    "GGT increased", 100,   0,      NA,    NA,       23,
    # AVAL missing cannot grade
    "GGT increased", NA,    0,      NA,    NA,       24,
  ) %>%
    # set BASE to be abnormal (HIGH HIGH) and create FLAG
    # set high_indicator to "HIGH HIGH" also in function call
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_ggt_ctcv5 <- exp_ggt_ctcv5_norm %>%
    bind_rows(exp_ggt_ctcv5_abn)
  test_high(expected = exp_ggt_ctcv5, meta = atoxgr_criteria_ctcv5, high = "HIGH HIGH")
})

### GGT increased (NCICTCv6)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >10.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >3.0 - 10.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR 1.0 - 1.5 x BL if BL was abnormal

## Test 30b: CTCAEv6 GGT increased ----
test_that("derive_var_atoxgr_dir Test 30b: CTCAEv6 GGT increased", {
  local_atoxgr_test_data()
  # V5 and V4 criteria identical when BASELINE normal
  exp_ggt_ctcv6_norm <- exp_ggt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  exp_ggt_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRLO, ~BASE, ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    0,      40,    NA,       13,
    NA_character_,   60,    0,      40,    NA,       14,
    "GGT increased", 401,   0,      40,    "4",      15,
    "GGT increased", 400,   0,      40,    "3",      16,
    "GGT increased", 121,   0,      40,    "3",      17,
    "GGT increased", 120,   0,      40,    "2",      18,
    "GGT increased", 61,    0,      40,    "2",      19,
    "GGT increased", 60,    0,      40,    "1",      20,
    "GGT increased", 40,    0,      40,    "1",      21,
    "GGT increased", 39,    0,      40,    "0",      22,
    # ANRHI missing - cannot grade
    "GGT increased", 100,   0,      NA,    NA,       23,
    # AVAL missing cannot grade
    "GGT increased", NA,    0,      NA,    NA,       24,
  ) %>%
    # set BASE to be abnormal (HIGH HIGH) and create FLAG
    # set high_indicator to "HIGH HIGH" also in function call
    mutate(
      AVALU = NA_character_,
      ANRHI = BASE - 1,
      BNRIND = "HIGH HIGH"
    )

  # combine records with baseline normal and abnormal
  exp_ggt_ctcv6 <- exp_ggt_ctcv6_norm %>%
    bind_rows(exp_ggt_ctcv6_abn)

  test_high(expected = exp_ggt_ctcv6, meta = atoxgr_criteria_ctcv6, high = "HIGH HIGH")
})

### Haptoglobin decreased (NCICTCv4)
# Same as NCICTCv5
### Grade 1: <LLN


## Test 31: CTCAEv4 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 31: CTCAEv4 Haptoglobin decreased", {
  local_atoxgr_test_data()
  test_low(expected = exp_hapt, meta = atoxgr_criteria_ctcv4)
})

## Test 32a: CTCAEv5 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 32a: CTCAEv5 Haptoglobin decreased", {
  local_atoxgr_test_data()
  test_low(expected = exp_hapt, meta = atoxgr_criteria_ctcv5)
})

## Test 32b: CTCAEv6 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 32b: CTCAEv6 Haptoglobin decreased", {
  local_atoxgr_test_data()
  test_low(expected = exp_hapt, meta = atoxgr_criteria_ctcv6)
})

### Hemoglobin increased
# NCICTCAEv5 same as NCICTCAEv4 when BASE is normal
# SI unit is "g/L" USCV unit is "g/dL"
### Grade 3: Increase in >4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 2: Increase in >2 - 4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 1: Increase in >0 - 2 gm/dL above ULN or above baseline if baseline is above ULN



## Test 33a: CTCAEv4 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 33a: CTCAEv4 Hemoglobin increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_hgbi_si, meta = atoxgr_criteria_ctcv4)
})

## Test 33b: CTCAEv4 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 33b: CTCAEv4 Hemoglobin increased (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_hgbi_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 34a: CTCAEv5 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 34a: CTCAEv5 Hemoglobin increased (SI unit)", {
  local_atoxgr_test_data()
  exp_hgbi_si <- exp_hgbi_si %>%
    filter(V5 == "Y")

  test_high(expected = exp_hgbi_si, meta = atoxgr_criteria_ctcv5)
})

## Test 34b: CTCAEv5 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 34b: CTCAEv5 Hemoglobin increased (USCV unit)", {
  local_atoxgr_test_data()
  exp_hgbi_cv <- exp_hgbi_cv %>%
    filter(V5 == "Y")

  test_high(expected = exp_hgbi_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 34c: CTCAEv6 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 34c: CTCAEv6 Hemoglobin increased (SI unit)", {
  local_atoxgr_test_data()
  exp_hgbi_si <- exp_hgbi_si %>%
    filter(V5 == "Y")

  test_high(expected = exp_hgbi_si, meta = atoxgr_criteria_ctcv6)
})

## Test 34d: CTCAEv6 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 34d: CTCAEv6 Hemoglobin increased (USCV unit)", {
  local_atoxgr_test_data()
  exp_hgbi_cv <- exp_hgbi_cv %>%
    filter(V5 == "Y")

  test_high(expected = exp_hgbi_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### INR increased (NCICTCV4)
### NCICTCV5 different for grade 1
### Grade 3: >2.5 x ULN; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5 x ULN; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1 - 1.5 x ULN; >1 - 1.5 times above baseline if on anticoagulation

## Test 35: CTCAEv4 INR increased ----
test_that("derive_var_atoxgr_dir Test 35: CTCAEv4 INR increased", {
  local_atoxgr_test_data()
  exp_inri <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~BASE, ~ANRHI, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    120,   200,    NA_character_, NA,       1,
    NA_character_,   60,    50,    100,    NA_character_, NA,       2,
    # GRADE derived from AVAL against ANRHI
    "INR IncreaSed", 251,   200,   100,    NA_character_, "3",      3,
    "INR Increased", 250,   199,   100,    NA_character_, "2",      4,
    "INR Increased", 151,   150,   100,    NA_character_, "2",      5,
    "INR Increased", 150,   150,   100,    NA_character_, "1",      6,
    "INR Increased", 101,   150,   100,    NA_character_, "1",      7,
    "INR Increased", 100,   100,   100,    NA_character_, "0",      8,
    # GRADE derived from AVAL against BASE
    "INR IncreaSed", 251,   100,   200,    NA_character_, "3",      9,
    "INR Increased", 250,   100,   199,    NA_character_, "2",      10,
    "INR Increased", 151,   100,   150,    NA_character_, "2",      11,
    "INR Increased", 150,   100,   150,    NA_character_, "1",      12,
    "INR Increased", 101,   100,   150,    NA_character_, "1",      13,
    "INR Increased", 100,   100,   100,    NA_character_, "0",      14,
    # BASE missing - AVAL <= ANRLO cannot grade as NORMAL
    "INR Increased", 100,   NA,    100,    NA_character_, NA,       15,
    # ANRHI missing - AVAL <= BASE cannot grade as NORMAL
    "INR Increased", 100,   100,   NA,     NA_character_, NA,       16,
    # AVAL missing cannot grade
    "INR Increased", NA,    100,   100,    NA_character_, NA,       17,
  )

  test_high(expected = exp_inri, meta = atoxgr_criteria_ctcv4)
})

### INR increased (NCICTCV5)
### NCICTCV5 different to NCICTCv4 (do not use x ULN)
### Grade 3: >2.5; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1.2 - 1.5; >1 - 1.5 times above baseline if on anticoagulation

## Test 36a: CTCAEv5 INR increased ----
test_that("derive_var_atoxgr_dir Test 36a: CTCAEv5 INR increased", {
  local_atoxgr_test_data()
  exp_inri <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~BASE, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    120,   NA_character_, NA,       1,
    NA_character_,   60,    50,    NA_character_, NA,       2,
    # GRADE derived from AVAL against first half of criteria
    "INR IncreaSed", 2.51,  2.6,   NA_character_, "3",      3,
    "INR Increased", 2.5,   1,     NA_character_, "2",      4,
    "INR Increased", 1.51,  1,     NA_character_, "2",      5,
    "INR Increased", 1.5,   1,     NA_character_, "1",      6,
    "INR Increased", 1.2,   1,     NA_character_, "1",      7,
    "INR Increased", 1.19,  1.19,  NA_character_, "0",      8,
    # GRADE derived from AVAL against BASE
    "INR IncreaSed", 2.5,   0.99,  NA_character_, "3",      9,
    "INR Increased", 1.5,   0.6,   NA_character_, "2",      10,
    "INR Increased", 1.5,   0.99,  NA_character_, "2",      11,
    "INR Increased", 1.2,   0.8,   NA_character_, "1",      12,
    "INR Increased", 1.2,   1.19,  NA_character_, "1",      13,
    "INR Increased", 1.2,   1.2,   NA_character_, "0",      14,
    # BASE missing - AVAL <= 1.2 cannot grade as NORMAL
    "INR Increased", 1.2,   NA,    NA_character_, NA,       15,
    # AVAL missing cannot grade
    "INR Increased", NA,    100,   NA_character_, NA,       16,
  )

  test_high(expected = exp_inri, meta = atoxgr_criteria_ctcv5)
})

## Test 36b: CTCAEv6 INR increased ----
test_that("derive_var_atoxgr_dir Test 36b: CTCAEv6 INR increased", {
  local_atoxgr_test_data()
  exp_inri <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~BASE, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    120,   NA_character_, NA,       1,
    NA_character_,   60,    50,    NA_character_, NA,       2,
    # GRADE derived from AVAL against first half of criteria
    "INR IncreaSed", 2.51,  2.6,   NA_character_, "3",      3,
    "INR Increased", 2.5,   1,     NA_character_, "2",      4,
    "INR Increased", 1.51,  1,     NA_character_, "2",      5,
    "INR Increased", 1.5,   1,     NA_character_, "1",      6,
    "INR Increased", 1.2,   1,     NA_character_, "1",      7,
    "INR Increased", 1.19,  1.19,  NA_character_, "0",      8,
    # GRADE derived from AVAL against BASE
    "INR IncreaSed", 2.5,   0.99,  NA_character_, "3",      9,
    "INR Increased", 1.5,   0.6,   NA_character_, "2",      10,
    "INR Increased", 1.5,   0.99,  NA_character_, "2",      11,
    "INR Increased", 1.2,   0.8,   NA_character_, "1",      12,
    "INR Increased", 1.2,   1.19,  NA_character_, "1",      13,
    "INR Increased", 1.2,   1.2,   NA_character_, "0",      14,
    # BASE missing - AVAL <= 1.2 cannot grade as NORMAL
    "INR Increased", 1.2,   NA,    NA_character_, NA,       15,
    # AVAL missing cannot grade
    "INR Increased", NA,    100,   NA_character_, NA,       16,
  )

  test_high(expected = exp_inri, meta = atoxgr_criteria_ctcv6)
})

### Lipase increased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN


## Test 37: CTCAEv4 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 37: CTCAEv4 Lipase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_lip, meta = atoxgr_criteria_ctcv4)
})

## Test 38a: CTCAEv5 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 38a: CTCAEv5 Lipase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_lip, meta = atoxgr_criteria_ctcv5)
})

### Lipase increased NCICTCAEv6 (small difference to v4 and v5)
### Grade 4: >5.0 x ULN
### Grade 3: >3.0 - 5.0 x ULN
### Grade 2: >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN


## Test 38b: CTCAEv6 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 38b: CTCAEv6 Lipase increased", {
  local_atoxgr_test_data()
  test_high(expected = exp_lipv6, meta = atoxgr_criteria_ctcv6)
})

### Lymphocyte count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <0.2 x 10e9 /L
### Grade 3: <0.5 - 0.2 x 10e9 /L
### Grade 2: <0.8 - 0.5 x 10e9 /L
### Grade 1: <LLN - 0.8 x 10e9/L


## Test 39a: CTCAEv4 Lymphocyte count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 39a: CTCAEv4 Lymphocyte count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_si, meta = atoxgr_criteria_ctcv4)
})


## Test 39b: CTCAEv4 Lymphocyte count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 39b: CTCAEv4 Lymphocyte count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 39c: CTCAEv4 Lymphocyte count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 39c: CTCAEv4 Lymphocyte count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_cv2, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 40a: CTCAEv5 Lymphocyte count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 40a: CTCAEv5 Lymphocyte count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_si, meta = atoxgr_criteria_ctcv5)
})

## Test 40b: CTCAEv5 Lymphocyte count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 40b: CTCAEv5 Lymphocyte count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_cv, meta = atoxgr_criteria_ctcv5_uscv)
})


## Test 40c: CTCAEv5 Lymphocyte count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 40c: CTCAEv5 Lymphocyte count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymd_cv2, meta = atoxgr_criteria_ctcv5_uscv)
})

### Lymphocyte count increased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 3: >20,000/mm3
### Grade 2: >4000/mm3 - 20,000/mm3


## Test 41a: CTCAEv4 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 41a: CTCAEv4 Lymphocyte count increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_si, meta = atoxgr_criteria_ctcv4)
})




## Test 41b: CTCAEv4 Lymphocyte count increased  (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 41b: CTCAEv4 Lymphocyte count increased  (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 42a: CTCAEv5 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 42a: CTCAEv5 Lymphocyte count increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_si, meta = atoxgr_criteria_ctcv5)
})

## Test 42b: CTCAEv5 Lymphocyte count increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 42b: CTCAEv5 Lymphocyte count increased (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 42c: CTCAEv6 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 42c: CTCAEv6 Lymphocyte count increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_si, meta = atoxgr_criteria_ctcv6)
})

## Test 42d: CTCAEv6 Lymphocyte count increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 42d: CTCAEv6 Lymphocyte count increased (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_lymi_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Neutrophil count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <1.0 - 0.5 x 10e9 /L
### Grade 2: <1.5 - 1.0 x 10e9 /L
### Grade 1: <LLN - 1.5 x 10e9 /L


## Test 43a: CTCAEv4 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 43a: CTCAEv4 Neutrophil count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_neut_si, meta = atoxgr_criteria_ctcv4)
})




## Test 43b: CTCAEv4 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 43b: CTCAEv4 Neutrophil count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_neut_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 44a: CTCAEv5 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 44a: CTCAEv5 Neutrophil count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_neut_si, meta = atoxgr_criteria_ctcv5)
})

## Test 44b: CTCAEv5 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44b: CTCAEv5 Neutrophil count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_neut_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

### Neutrophil count decreased NCICTCAEv6 (different to v4 and v5)
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <0.1 x 10^9/L
### Grade 3: <0.5 - 0.1 x 10^9/L
### Grade 2: <1.0 - 0.5 x 10^9/L
### Grade 1: <1.5 - 1.0 x 10^9/L


## Test 44c: CTCAEv6 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 44c: CTCAEv6 Neutrophil count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_neut_siv6, meta = atoxgr_criteria_ctcv6)
})

## Test 44d: CTCAEv6 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44d: CTCAEv6 Neutrophil count decreased (USCV unit)", {
  local_atoxgr_test_data()
  exp_neut_usv6 <- exp_neut_siv6 %>%
    mutate(
      AVALU = if_else(AVALU == "10^9/L", "10^3/uL", AVALU)
    )

  test_low(expected = exp_neut_usv6, meta = atoxgr_criteria_ctcv6_uscv)
})

## Test 44e: CTCAEv6 Neutrophil count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44e: CTCAEv6 Neutrophil count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  exp_neut_us2v6 <- exp_neut_siv6 %>%
    mutate(
      AVALU = if_else(AVALU == "10^9/L", "10^3/mL", AVALU),
      AVAL = AVAL * 1000
    )

  test_low(expected = exp_neut_us2v6, meta = atoxgr_criteria_ctcv6_uscv)
})

### Platelet count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <50.0 - 25.0 x 10e9 /L
### Grade 2: <75.0 - 50.0 x 10e9 /L
### Grade 1: <LLN - 75.0 x 10e9 /L


## Test 45a: CTCAEv4 Platelet count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 45a: CTCAEv4 Platelet count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_si, meta = atoxgr_criteria_ctcv4)
})


## Test 45b: CTCAEv4 Platelet count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 45b: CTCAEv4 Platelet count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 45c: CTCAEv4 Platelet count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 45c: CTCAEv4 Platelet count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_cv2, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 46a: CTCAEv5 Platelet count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 46a: CTCAEv5 Platelet count decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_si, meta = atoxgr_criteria_ctcv5)
})

## Test 46b: CTCAEv5 Platelet count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46b: CTCAEv5 Platelet count decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 46c: CTCAEv5 Platelet count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46c: CTCAEv5 Platelet count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_cv2, meta = atoxgr_criteria_ctcv5_uscv)
})

### For NCICTCAEv6 Term: Thrombocytopenia replaces Term: Platelet count decreased in NCICTCAEv5
### NCICTCAEv5 and NCICTCAEv6 criteria is the same apart from grade 4
### For Grade 4: it was <25.0 x 10e9 /L but is now <10.0 x 10e9 /L
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <10.0 x 10e9 /L
### Grade 3: <50.0 - 25.0 x 10e9 /L
### Grade 2: <75.0 - 50.0 x 10e9 /L
### Grade 1: <LLN - 75.0 x 10e9 /L


## Test 46d: CTCAEv6 Thrombocytopenia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 46d: CTCAEv6 Thrombocytopenia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plate_v6_si, meta = atoxgr_criteria_ctcv6)
})

## Test 46e: CTCAEv6 Thrombocytopenia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46e: CTCAEv6 Thrombocytopenia (USCV unit)", {
  local_atoxgr_test_data()
  exp_plate_v6_cv <- exp_plate_v6_si %>%
    mutate(AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU))

  test_low(expected = exp_plate_v6_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

## Test 46f: CTCAEv6 Thrombocytopenia (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46f: CTCAEv6 Platelet count decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  exp_plate_v6cv2 <- exp_plate_v6_si %>%
    mutate(
      AVAL = AVAL * 1000,
      ANRLO = ANRLO * 1000,
      ANRHI = ANRHI * 1000,
      AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
    )

  test_low(expected = exp_plate_v6cv2, meta = atoxgr_criteria_ctcv6_uscv)
})

### Serum amylase increased
### NCICTCAEv4 and NCICTCAEv5 criteria essentially the same
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN


## Test 47a: CTCAEv4 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 47a: CTCAEv4 Serum amylase increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_seri, meta = atoxgr_criteria_ctcv4)
})

## Test 47b: CTCAEv4 Serum amylase increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 47b: CTCAEv4 Serum amylase increased (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_seri, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 48a: CTCAEv5 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 48a: CTCAEv5 Serum amylase Nincreased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_seri, meta = atoxgr_criteria_ctcv5)
})

## Test 48b: CTCAEv5 Serum amylase increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 48b: CTCAEv5 Serum amylase increased (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_seri, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 48c: CTCAEv6 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 48c: CTCAEv6 Serum amylase increased (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_seri, meta = atoxgr_criteria_ctcv6)
})

### White blood cell decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <1.0 x 10e9 /L
### Grade 3: <2.0 - 1.0 x 10e9 /L
### Grade 2: <3.0 - 2.0 x 10e9 /L
### Grade 1: <LLN - 3.0 x 10e9 /L





## Test 49a: CTCAEv4 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 49a: CTCAEv4 White blood cell decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd, meta = atoxgr_criteria_ctcv4)
})

## Test 49b: CTCAEv4 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 49b: CTCAEv4 White blood cell decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_uscv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 50a: CTCAEv5 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 50a: CTCAEv5 White blood cell decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd, meta = atoxgr_criteria_ctcv5)
})

## Test 50b: CTCAEv5 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 50b: CTCAEv5 White blood cell decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_uscv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 50c: CTCAEv6 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 50c: CTCAEv6 White blood cell decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd, meta = atoxgr_criteria_ctcv6)
})

## Test 50d: CTCAEv6 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 50d: CTCAEv6 White blood cell decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_uscv, meta = atoxgr_criteria_ctcv6_uscv)
})

## Metabolism and nutrition disorders

### Hypercalcemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit
### Grade 4: >3.4 mmol/L
### Grade 3: >3.1 - 3.4 mmol/L
### Grade 2: >2.9 - 3.1 mmol/L
### Grade 1: >ULN - 2.9 mmol/L

### CV unit
### Grade 4: >13.5 mg/dL
### Grade 3: >12.5 - 13.5 mg/dL
### Grade 2: >11.5 - 12.5 mg/dL
### Grade 1: >ULN - 11.5 mg/dL



## Test 51a: CTCAEv4 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 51a: CTCAEv4 Hypercalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_si, meta = atoxgr_criteria_ctcv4)
})

## Test 51b: CTCAEv4 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 51b: CTCAEv4 Hypercalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 52a: CTCAEv5 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 52a: CTCAEv5 Hypercalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_si, meta = atoxgr_criteria_ctcv5)
})

## Test 52b: CTCAEv5 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 52b: CTCAEv5 Hypercalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 52c: CTCAEv6 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 52c: CTCAEv6 Hypercalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_si, meta = atoxgr_criteria_ctcv6)
})

## Test 52c: CTCAEv6 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 52c: CTCAEv6 Hypercalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypercalcemia (Ionized)
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit
### Grade 4: >1.8 mmol/L
### Grade 3: >1.6 - 1.8 mmol/L
### Grade 2: >1.5 - 1.6 mmol/L
### Grade 1: >ULN - 1.5 mmol/L

### CV unit
### Grade 4: >7.2 mg/dL
### Grade 3: >6.4 - 7.2 mg/dL
### Grade 2: >6 - 6.4 mg/dL
### Grade 1: >ULN - 11.5 mg/dL



## Test 53a: CTCAEv4 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 53a: CTCAEv4 Hypercalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_si, meta = atoxgr_criteria_ctcv4)
})

## Test 53b: CTCAEv4 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 53b: CTCAEv4 Hypercalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 54a: CTCAEv5 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 54a: CTCAEv5 Hypercalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_si, meta = atoxgr_criteria_ctcv5)
})

## Test 54b: CTCAEv5 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 54b: CTCAEv5 Hypercalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 54c: CTCAEv6 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 54c: CTCAEv6 Hypercalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_si, meta = atoxgr_criteria_ctcv6)
})

## Test 54d: CTCAEv6 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 54d: CTCAEv6 Hypercalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hyperglycemia (Fasting) (NCICTCAEv4 and NCICTCAEv6)
### not included in NCICTCAEv5
### SI unit
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L
### Grade 2: >8.9 - 13.9 mmol/L
### Grade 1: >ULN - 8.9 mmol/L

### CV unit
### Grade 4: >500 mg/dL
### Grade 3: >250 - 500 mg/dL
### Grade 2: >160 - 250 mg/dL
### Grade 1: >ULN - 160 mg/dL


## Test 55a: CTCAEv4 Hyperglycemia (Fasting) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 55a: CTCAEv4 Hyperglycemia (Fasting) (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glycfi, meta = atoxgr_criteria_ctcv4)
})

## Test 55b: CTCAEv4 Hyperglycemia (Fasting) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 55b: CTCAEv4 Hyperglycemia (Fasting) (USCV unit)", {
  local_atoxgr_test_data()
  exp_glycfi <- tibble::tribble(
    ~ATOXDSCH,                  ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU, ~ATOXGRH,  ~TESTNUM,
    "Not a term",               501,    0,       100,     "mg/dL",       NA,  1,
    NA_character_,              501,    0,       100,     "mg/dL",       NA,  2,
    # ANRHI not missing
    "Hyperglycemia (Fasting)",  501,    0,       100,     "mg/dL",      "4",  3,
    "Hyperglycemia (Fasting)",  500,    0,       100,     "mg/dL",      "3",  4,
    "Hyperglycemia (Fasting)",  251,    0,       100,     "mg/dL",      "3",  5,
    "Hyperglycemia (Fasting)",  250,    0,       100,     "mg/dL",      "2",  6,
    "Hyperglycemia (Fasting)",  161,    0,       100,     "mg/dL",      "2",  7,
    "Hyperglycemia (Fasting)",  160,    0,       100,     "mg/dL",      "1",  8,
    "Hyperglycemia (Fasting)",  101,    0,       100,     "mg/dL",      "1",  9,
    "Hyperglycemia (Fasting)",  100,    0,       100,     "mg/dL",      "0",  10,
    # ANRHI missing - can grade 2-4
    "Hyperglycemia (Fasting)",  501,    0,       NA,      "mg/dL",      "4",  11,
    "Hyperglycemia (Fasting)",  500,    0,       NA,      "mg/dL",      "3",  12,
    "Hyperglycemia (Fasting)",  251,    0,       NA,      "mg/dL",      "3",  13,
    "Hyperglycemia (Fasting)",  250,    0,       NA,      "mg/dL",      "2",  14,
    "Hyperglycemia (Fasting)",  161,    0,       NA,      "mg/dL",      "2",  15,
    # ANRHI missing - can NOT grade 0 or 1
    "Hyperglycemia (Fasting)",  160,    0,       NA,      "mg/dL",       NA,  16,
    "Hyperglycemia (Fasting)",  101,    0,       NA,      "mg/dL",       NA,  17,
    "Hyperglycemia (Fasting)",  100,    0,       NA,      "mg/dL",       NA,  18,
    # Unit missing cannot grade
    "Hyperglycemia (Fasting)",  100,    0,       5.3,          NA,       NA,  19,
    # AVAL missing cannot grade
    "Hyperglycemia (Fasting)",  NA,     0,       5.3,     "mg/dL",       NA,  20,
  )

  test_high(expected = exp_glycfi, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 55c: CTCAEv6 Hyperglycemia (Fasting) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 55c: CTCAEv6 Hyperglycemia (Fasting) (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glycfi, meta = atoxgr_criteria_ctcv6)
})

### Hyperglycemia (NCICTCAEv4 and NCICTCAEv6)
### not included in NCICTCAEv5
### SI unit
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L

### CV unit
### Grade 4: >500 mg/dL
### Grade 3: >250 - 500 mg/dL



## Test 56a: CTCAEv4 Hyperglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 56a: CTCAEv4 Hyperglycemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glyci, meta = atoxgr_criteria_ctcv4)
})

## Test 56b: CTCAEv4 Hyperglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 56b: CTCAEv4 Hyperglycemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glyci_uscv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 56c: CTCAEv6 Hyperglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 56c: CTCAEv6 Hyperglycemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glyci, meta = atoxgr_criteria_ctcv6)
})

## Test 56d: CTCAEv6 Hyperglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 56d: CTCAEv6 Hyperglycemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_glyci_uscv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hyperkalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI + CV unit
### Grade 4: >7.0 mmol/L
### Grade 3: >6.0 - 7.0 mmol/L
### Grade 2: >5.5 - 6.0 mmol/L
### Grade 1: >ULN - 5.5 mmol/L


## Test 57a: CTCAEv4 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 57a: CTCAEv4 Hyperkalemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_kalei, meta = atoxgr_criteria_ctcv4)
})

## Test 57b: CTCAEv4 Hyperkalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 57b: CTCAEv4 Hyperkalemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_kalei, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 58a: CTCAEv5 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 58a: CTCAEv5 Hyperkalemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_kalei, meta = atoxgr_criteria_ctcv5)
})

## Test 58b: CTCAEv5 Hyperkalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 58b: CTCAEv5 Hyperkalemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_kalei, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 58c: CTCAEv6 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 58c: CTCAEv6 Hyperkalemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_kalei, meta = atoxgr_criteria_ctcv6)
})

### Hypermagnesemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit
### Grade 4: >3.30 mmol/L
### Grade 3: >1.23 - 3.30 mmol/L
### Grade 1: >ULN - 1.23 mmol/L


### CV unit
### Grade 4: >8 mg/dL
### Grade 3: >3 - 8 mg/dL
### Grade 1: >ULN - 3 mg/dL



## Test 59a: CTCAEv4 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 59a: CTCAEv4 Hypermagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_si, meta = atoxgr_criteria_ctcv4)
})

## Test 59b: CTCAEv4 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 59b: CTCAEv4 Hypermagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 60a: CTCAEv5 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 60a: CTCAEv5 Hypermagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_si, meta = atoxgr_criteria_ctcv5)
})

## Test 60b: CTCAEv5 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 60b: CTCAEv5 Hypermagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 60c: CTCAEv6 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 60c: CTCAEv6 Hypermagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_si, meta = atoxgr_criteria_ctcv6)
})

## Test 60d: CTCAEv6 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 60d: CTCAEv6 Hypermagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_magni_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypernatremia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI + CV unit
### Grade 4: >160 mmol/L
### Grade 3: >155 - 160 mmol/L
### Grade 2: >150 - 155 mmol/L
### Grade 1: >ULN - 150 mmol/L


## Test 61a: CTCAEv4 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 61a: CTCAEv4 Hypernatremia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_natri, meta = atoxgr_criteria_ctcv4)
})

## Test 61b: CTCAEv4 Hypernatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 61b: CTCAEv4 Hypernatremia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_natri, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 62a: CTCAEv5 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 62a: CTCAEv5 Hypernatremia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_natri, meta = atoxgr_criteria_ctcv5)
})

## Test 62b: CTCAEv5 Hypernatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 62b: CTCAEv5 Hypernatremia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_natri, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 62c: CTCAEv6 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 62c: CTCAEv6 Hypernatremia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_natri, meta = atoxgr_criteria_ctcv6)
})

### Hypertriglyceridemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit
### Grade 4: >11.4 mmol/L
### Grade 3: >5.7 mmol/L - 11.4 mmol/L
### Grade 2: >3.42 mmol/L - 5.7 mmol/L
### Grade 1: 1.71 mmol/L - 3.42 mmol/L


### CV unit
### Grade 4: >1000 mg/dL
### Grade 3: >500 - 1000 mg/dL
### Grade 2: >300 - 500 mg/dL
### Grade 1: 150 - 300 mg/dL



## Test 63a: CTCAEv4 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 63a: CTCAEv4 Hypertriglyceridemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_si, meta = atoxgr_criteria_ctcv4)
})

## Test 63b: CTCAEv4 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 63b: CTCAEv4 Hypertriglyceridemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 64a: CTCAEv5 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 64a: CTCAEv5 Hypertriglyceridemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_si, meta = atoxgr_criteria_ctcv5)
})

## Test 64b: CTCAEv5 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 64b: CTCAEv5 Hypertriglyceridemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 64c: CTCAEv6 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 64c: CTCAEv6 Hypertriglyceridemia (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_si, meta = atoxgr_criteria_ctcv6)
})

## Test 64d: CTCAEv6 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 64d: CTCAEv6 Hypertriglyceridemia (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_trigi_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hyperuricemia (NCICTCAEv4)
### NCICTCAEv5 only has grade 3
### SI unit is umol/L = 1000 * mmol/L
### Grade 4: >0.59 mmol/L;
### Grade 3: >ULN - 0.59 mmol/L

### CV unit is mg/dL
### Grade 4: >10 mg/dL
### Grade 3: >ULN - 10 mg/dL




### Hyperuricemia (NCICTCAEv5)
### NCICTCAEv5 only has grade 3
### Grade 3: >ULN

## Test 65a: CTCAEv5 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 65a: CTCAEv5 Hyperuricemia (SI unit)", {
  local_atoxgr_test_data()
  exp_urici_si <- exp_urici_si %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")

  test_high(expected = exp_urici_si, meta = atoxgr_criteria_ctcv5)
})

## Test 65b: CTCAEv5 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 65b: CTCAEv5 Hyperuricemia (USCV unit)", {
  local_atoxgr_test_data()
  exp_urici_cv <- exp_urici_cv %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")

  test_high(expected = exp_urici_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 65c: CTCAEv6 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 65c: CTCAEv6 Hyperuricemia (SI unit)", {
  local_atoxgr_test_data()
  exp_urici_si <- exp_urici_si %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")

  test_high(expected = exp_urici_si, meta = atoxgr_criteria_ctcv6)
})

## Test 65d: CTCAEv6 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 65d: CTCAEv6 Hyperuricemia (USCV unit)", {
  local_atoxgr_test_data()
  exp_urici_cv <- exp_urici_cv %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")

  test_high(expected = exp_urici_cv, meta = atoxgr_criteria_ctcv6_uscv)
})


# If unit missing then grade CANNOT be calculated as needed for grade 4
## Test 66a: CTCAEv4 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 66a: CTCAEv4 Hyperuricemia (SI unit)", {
  local_atoxgr_test_data()
  exp_urici_si <- exp_urici_si %>%
    mutate(ATOXGRH = if_else(is.na(AVALU), NA_character_, ATOXGRH))

  input_urici_si <- exp_urici_si %>%
    select(-ATOXGRH)

  test_high(expected = exp_urici_si, meta = atoxgr_criteria_ctcv4)
})

## Test 66b: CTCAEv4 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 66b: CTCAEv4 Hyperuricemia (USCV unit)", {
  local_atoxgr_test_data()
  exp_urici_cv <- exp_urici_cv %>%
    mutate(ATOXGRH = if_else(is.na(AVALU), NA_character_, ATOXGRH))

  test_high(expected = exp_urici_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

### Hypoalbuminemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI init is g/L - CV unit is g/dL
### Grade 3: <20 g/L
### Grade 2: <30 - 20 g/L
### Grade 1: <LLN - 30 g/L



## Test 67a: CTCAEv4 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 67a: CTCAEv4 Hypoalbuminemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_si, meta = atoxgr_criteria_ctcv4)
})

## Test 67b: CTCAEv4 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 67b: CTCAEv4 Hypoalbuminemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 68a: CTCAEv5 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 68a: CTCAEv5 Hypoalbuminemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_si, meta = atoxgr_criteria_ctcv5)
})

## Test 68b: CTCAEv5 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 68b: CTCAEv5 Hypoalbuminemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 68c: CTCAEv6 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 68c: CTCAEv6 Hypoalbuminemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_si, meta = atoxgr_criteria_ctcv6)
})

## Test 68d: CTCAEv6 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 68d: CTCAEv6 Hypoalbuminemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albd_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypocalcemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is mmol/L
### Grade 4: <1.5 mmol/L
### Grade 3: <1.75 - 1.5 mmol/L
### Grade 2: <2.0 - 1.75 mmol/L
### Grade 1: <LLN - 2.0 mmol/L

### CV unit is mg/dL
### Grade 4: <6 mg/dL
### Grade 3: <7 - 6 mg/dL
### Grade 2: <8 - 7 mg/dL
### Grade 1: <LLN - 8 mg/dL



## Test 69a: CTCAEv4 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 69a: CTCAEv4 Hypocalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_si, meta = atoxgr_criteria_ctcv4)
})

## Test 69b: CTCAEv4 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 69b: CTCAEv4 Hypocalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 70a: CTCAEv5 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 70a: CTCAEv5 Hypocalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_si, meta = atoxgr_criteria_ctcv5)
})

## Test 70b: CTCAEv5 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 70b: CTCAEv5 Hypocalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 70c: CTCAEv6 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 70c: CTCAEv6 Hypocalcemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_si, meta = atoxgr_criteria_ctcv6)
})

## Test 70d: CTCAEv6 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 70d: CTCAEv6 Hypocalcemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypocalcemia (Ionized)
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is mmol/L
### Grade 4: <0.8 mmol/L
### Grade 3: <0.9 - 0.8 mmol/L
### Grade 2: <1.0 - 0.9 mmol/L
### Grade 1: <LLN - 1.0 mmol/L

### CV unit is mg/dL
### mg/dL is 4 x mmol/L




## Test 71a: CTCAEv4 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 71a: CTCAEv4 Hypocalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_si, meta = atoxgr_criteria_ctcv4)
})

## Test 71b: CTCAEv4 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 71b: CTCAEv4 Hypocalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 72a: CTCAEv5 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 72a: CTCAEv5 Hypocalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_si, meta = atoxgr_criteria_ctcv5)
})

## Test 72b: CTCAEv5 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 72b: CTCAEv5 Hypocalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 72c: CTCAEv6 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 72c: CTCAEv6 Hypocalcemia (Ionized) (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_si, meta = atoxgr_criteria_ctcv6)
})

## Test 72d: CTCAEv6 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 72d: CTCAEv6 Hypocalcemia (Ionized) (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypoglycemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is mmol/L
### Grade 4: <1.7 mmol/L
### Grade 3: <2.2 - 1.7 mmol/L
### Grade 2: <3.0 - 2.2 mmol/L
### Grade 1: <LLN - 3.0 mmol/L


### CV unit is mg/dL
### Grade 4: <30 mg/dL
### Grade 3: <40 - 30 mg/dL
### Grade 2: <55 - 40 mg/dL
### Grade 1: <LLN - 55 mg/dL



## Test 73a: CTCAEv4 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 73a: CTCAEv4 Hypoglycemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_si, meta = atoxgr_criteria_ctcv4)
})

## Test 73b: CTCAEv4 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 73b: CTCAEv4 Hypoglycemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_cv, meta = atoxgr_criteria_ctcv4_uscv)
})


## Test 74a: CTCAEv5 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 74a: CTCAEv5 Hypoglycemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_si, meta = atoxgr_criteria_ctcv5)
})

## Test 74b: CTCAEv5 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 74b: CTCAEv5 Hypoglycemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 74c: CTCAEv6 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 74c: CTCAEv6 Hypoglycemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_si, meta = atoxgr_criteria_ctcv6)
})

## Test 74d: CTCAEv6 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 74d: CTCAEv6 Hypoglycemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glycd_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypokalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit and CV unit is mmol/L
### Grade 4: <2.5 mmol/L
### Grade 3: <3.0 - 2.5 mmol/L
### Grade 2: <LLN - 3.0 mmol/L


## Test 75a: CTCAEv4 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 75a: CTCAEv4 Hypokalemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv4)
})

## Test 75b: CTCAEv4 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 75b: CTCAEv4 Hypokalemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 76a: CTCAEv5 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 76a: CTCAEv5 Hypokalemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv5)
})

## Test 76b: CTCAEv5 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 76b: CTCAEv5 Hypokalemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 76c: CTCAEv6 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 76c: CTCAEv6 Hypokalemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv6)
})

## Test 76d: CTCAEv6 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 76d: CTCAEv6 Hypokalemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_kaled, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hypomagnesemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is mmol/L
### Grade 4: <0.3 mmol/L
### Grade 3: <0.4 - 0.3 mmol/L
### Grade 2: <0.5 - 0.4 mmol/L
### Grade 1: <LLN - 0.5 mmol/L


### CV unit is mg/dL
### Grade 4: <0.7 mg/dL
### Grade 3: <0.9 - 0.7 mg/dL
### Grade 2: <1.2 - 0.9 mg/dL
### Grade 1: <LLN - 1.2 mg/dL



## Test 77a: CTCAEv4 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 77a: CTCAEv4 Hypomagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_si, meta = atoxgr_criteria_ctcv4)
})

## Test 77b: CTCAEv4 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 77b: CTCAEv4 Hypomagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_cv, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 78a: CTCAEv5 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 78a: CTCAEv5 Hypomagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_si, meta = atoxgr_criteria_ctcv5)
})

## Test 78b: CTCAEv5 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 78b: CTCAEv5 Hypomagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_cv, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 78c: CTCAEv6 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 78c: CTCAEv6 Hypomagnesemia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_si, meta = atoxgr_criteria_ctcv6)
})

## Test 78d: CTCAEv6 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 78d: CTCAEv6 Hypomagnesemia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magnd_cv, meta = atoxgr_criteria_ctcv6_uscv)
})

### Hyponatremia (NCICTCAEv4)
### NCICTCAEv4 and NCICTCAEv5 essentially the same (slightly different text)
### SI unit and CV unit is mmol/L
### Grade 4: <120 mmol/L
### Grade 3: <130 - 120 mmol/L
### Grade 1: <LLN - 130 mmol/L


## Test 79a: CTCAEv4 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 79a: CTCAEv4 Hyponatremia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_natrd, meta = atoxgr_criteria_ctcv4)
})

## Test 79b: CTCAEv4 Hyponatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 79b: CTCAEv4 Hyponatremia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_natrd, meta = atoxgr_criteria_ctcv4_uscv)
})

## Test 80a: CTCAEv5 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 80a: CTCAEv5 Hyponatremia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_natrd, meta = atoxgr_criteria_ctcv5)
})

## Test 80b: CTCAEv5 Hyponatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 80b: CTCAEv5 Hyponatremia (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_natrd, meta = atoxgr_criteria_ctcv5_uscv)
})

## Test 80c: CTCAEv6 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 80c: CTCAEv6 Hyponatremia (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_natrd, meta = atoxgr_criteria_ctcv6)
})

### Hypophosphatemia
### Only in NCICTCAEv4 (not NCICTCAEv5)
### SI unit is mmol/L
### Grade 4: <0.3 mmol/L
### Grade 3: <0.6 - 0.3 mmol/L
### Grade 2: <0.8 - 0.6 mmol/L
### Grade 1: <LLN - 0.8 mmol/L


### CV unit is mg/dL
### Grade 4: <1 mg/dL
### Grade 3: <2 - 1 mg/dL
### Grade 2: <2.5 - 2 mg/dL
### Grade 1: <LLN - 2.5 mg/dL

## Test 81a: CTCAEv4 Hypophosphatemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 81: CTCAEv4 Hypophosphatemia (SI unit)", {
  local_atoxgr_test_data()
  exp_phosd <- tibble::tribble(
    ~ATOXDSCL,           ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
    "Not a term",        0.29,   1,      100,    "mmol/L",  NA,       1,
    NA_character_,       0.29,   1,      100,    "mmol/L",  NA,       2,
    # ANRLO not missing
    "Hypophosphatemia",  0.29,   1,      100,    "mmol/L",  "4",      3,
    "Hypophosphatemia",  0.3,    1,      100,    "mmol/L",  "3",      4,
    "Hypophosphatemia",  0.59,   1,      100,    "mmol/L",  "3",      5,
    "Hypophosphatemia",  0.6,    1,      100,    "mmol/L",  "2",      6,
    "Hypophosphatemia",  0.79,   1,      100,    "mmol/L",  "2",      7,
    "Hypophosphatemia",  0.8,    1,      100,    "mmol/L",  "1",      8,
    "Hypophosphatemia",  0.9,    1,      100,    "mmol/L",  "1",      9,
    "Hypophosphatemia",  1,      1,      100,    "mmol/L",  "0",      10,
    # ANRLO missing - can grade 3-4
    "Hypophosphatemia",  0.29,   NA,     100,    "mmol/L",  "4",      11,
    "Hypophosphatemia",  0.3,    NA,     100,    "mmol/L",  "3",      12,
    "Hypophosphatemia",  0.59,   NA,     100,    "mmol/L",  "3",      13,
    "Hypophosphatemia",  0.6,    NA,     100,    "mmol/L",  "2",      14,
    "Hypophosphatemia",  0.79,   NA,     100,    "mmol/L",  "2",      15,
    # ANRLO missing - can NOT grade 0 or 1
    "Hypophosphatemia",  0.8,    NA,     100,    "mmol/L",  NA,       16,
    "Hypophosphatemia",  0.9,    NA,     100,    "mmol/L",  NA,       17,
    "Hypophosphatemia",  1,      NA,     100,    "mmol/L",  NA,       18,
    # Unit missing cannot grade
    "Hypophosphatemia",  1,      1,      100,    NA,        NA,       19,
    # AVAL missing cannot grade
    "Hypophosphatemia",  NA,     1,      100,    "mmol/L",  NA,       20,
  )

  test_low(expected = exp_phosd, meta = atoxgr_criteria_ctcv4)
})

## Test 81b: CTCAEv4 Hypophosphatemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 81: CTCAEv4 Hypophosphatemia (USCV unit)", {
  local_atoxgr_test_data()
  exp_phosd <- tibble::tribble(
    ~ATOXDSCL,           ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL, ~TESTNUM,
    "Not a term",        0.99,   3,       100,     "mg/dL",        NA, 1,
    NA_character_,       0.99,   3,       100,     "mg/dL",        NA, 2,
    # ANRLO not missing
    "Hypophosphatemia",  0.99,   3,       100,     "mg/dL",       "4", 3,
    "Hypophosphatemia",  1,      3,       100,     "mg/dL",       "3", 4,
    "Hypophosphatemia",  1.9,    3,       100,     "mg/dL",       "3", 5,
    "Hypophosphatemia",  2,      3,       100,     "mg/dL",       "2", 6,
    "Hypophosphatemia",  2.49,   3,       100,     "mg/dL",       "2", 7,
    "Hypophosphatemia",  2.5,    3,       100,     "mg/dL",       "1", 8,
    "Hypophosphatemia",  2.9,    3,       100,     "mg/dL",       "1", 9,
    "Hypophosphatemia",  3,      3,       100,     "mg/dL",       "0", 10,
    # ANRLO missing - can grade 3-4
    "Hypophosphatemia",  0.99,   NA,      100,     "mg/dL",       "4", 11,
    "Hypophosphatemia",  1,      NA,      100,     "mg/dL",       "3", 12,
    "Hypophosphatemia",  1.9,    NA,      100,     "mg/dL",       "3", 13,
    "Hypophosphatemia",  2,      NA,      100,     "mg/dL",       "2", 14,
    "Hypophosphatemia",  2.49,   NA,      100,     "mg/dL",       "2", 15,
    # ANRLO missing - can NOT grade 0 or 1
    "Hypophosphatemia",  2.5,    NA,      100,     "mg/dL",        NA, 16,
    "Hypophosphatemia",  2.9,    NA,      100,     "mg/dL",        NA, 17,
    "Hypophosphatemia",  3,      NA,      100,     "mg/dL",        NA, 18,
    # Unit missing cannot grade
    "Hypophosphatemia",  3,      1,       100,          NA,        NA, 19,
    # AVAL missing cannot grade
    "Hypophosphatemia",  NA,     1,       100,     "mg/dL",        NA, 20,
  )

  test_low(expected = exp_phosd, meta = atoxgr_criteria_ctcv4_uscv)
})

# DAIDS grading----

### Acidosis
### Grade 4: pH < 7.3 with lifethreatening consequences
### Grade 3: pH < 7.3 without lifethreatening consequences
### Grade 2: pH >= 7.3 to < LLN

## Test 82: DAIDS Acidosis ----
test_that("derive_var_atoxgr_dir Test 82: DAIDS Acidosis", {
  local_atoxgr_test_data()
  exp_acido_daids <- tibble::tribble(
    ~ATOXDSCL,     ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRL, ~TESTNUM,
    "Not a term",  7.3,    7.35,   7.4,    NA,       1,
    NA_character_, 7.3,    7.35,   7.4,    NA,       2,
    # ANRLO not missing
    "Acidosis",    7.29,   7.35,   7.4,    "4",      3,
    "Acidosis",    7.3,    7.35,   7.4,    "2",      4,
    "Acidosis",    7.34,   7.35,   7.4,    "2",      5,
    "Acidosis",    7.35,   7.35,   7.4,    "0",      6,
    "Acidosis",    7.36,   7.35,   7.4,    "0",      7,
    # ANRLO missing - can grade 4
    "Acidosis",    7.29,   NA,     7.4,    "4",      8,
    # ANRLO missing - can NOT grade 0 or 2
    "Acidosis",    7.3,    NA,     7.4,    NA,       9,
    "Acidosis",    7.34,   NA,     7.4,    NA,       10,
    "Acidosis",    7.35,   NA,     7.4,    NA,       11,
    "Acidosis",    7.36,   NA,     7.4,    NA,       12,
    # AVAL missing cannot grade
    "Acidosis",    NA,     1.1,    1.4,    NA,       13,
  ) %>%
    mutate(AVALU = NA_character_)

  test_low(expected = exp_acido_daids, meta = atoxgr_criteria_daids)
})

### Albumin, Low
### SI unit is g/L CV unit is g/dL
### Grade 3: < 20
### Grade 2: >= 20 to < 30
### Grade 1: 30 to < LLN



## Test 83a: DAIDS Albumin, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 83a: DAIDS Albumin, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albl_daids_si, meta = atoxgr_criteria_daids)
})

## Test 83b: DAIDS Albumin, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 83b: DAIDS Albumin, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_albl_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Alkaline Phosphatase, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 84: DAIDS Alkaline Phosphatase, High ----
test_that("derive_var_atoxgr_dir Test 84: DAIDS Alkaline Phosphatase, High", {
  local_atoxgr_test_data()
  exp_alkp_i_daids <- tibble::tribble(
    ~ATOXDSCH,                    ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",                 30,     40,     NA,       1,
    NA_character_,                30,     40,     NA,       2,
    # ANRHI not missing
    "Alkaline Phosphatase, High", 401,    40,     "4",      3,
    "Alkaline Phosphatase, High", 400,    40,     "4",      4,
    "Alkaline Phosphatase, High", 399,    40,     "3",      5,
    "Alkaline Phosphatase, High", 200,    40,     "3",      6,
    "Alkaline Phosphatase, High", 199,    40,     "2",      7,
    "Alkaline Phosphatase, High", 100,    40,     "2",      8,
    "Alkaline Phosphatase, High", 99,     40,     "1",      9,
    "Alkaline Phosphatase, High", 51,     40,     "1",      10,
    "Alkaline Phosphatase, High", 50,     40,     "1",      11,
    "Alkaline Phosphatase, High", 49,     40,     "0",      12,
    # ANRHI missing cannot grade
    "Alkaline Phosphatase, High", 49,     NA,     NA,       13,
    # AVAL missing cannot grade
    "Alkaline Phosphatase, High", NA,     40,     NA,       14,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_alkp_i_daids, meta = atoxgr_criteria_daids)
})

### Alkalosis
### Grade 4: pH > 7.5 with lifethreatening consequences
### Grade 3: pH > 7.5 without lifethreatening consequences
### Grade 2: pH > ULN to <= 7.5

## Test 85: DAIDS Alkalosis ----
test_that("derive_var_atoxgr_dir Test 85: DAIDS Alkalosis", {
  local_atoxgr_test_data()
  exp_alkalo_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",  7.3,    7.35,   7.4,    NA,       1,
    NA_character_, 7.3,    7.35,   7.4,    NA,       2,
    # ANRHI not missing
    "Alkalosis",   7.51,   7.35,   7.4,    "4",      3,
    "Alkalosis",   7.5,    7.35,   7.4,    "2",      4,
    "Alkalosis",   7.41,   7.35,   7.4,    "2",      5,
    "Alkalosis",   7.4,    7.35,   7.4,    "0",      6,
    "Alkalosis",   7.39,   7.35,   7.4,    "0",      7,
    # ANRHI missing - can grade 4
    "Alkalosis",   7.51,   7.35,   NA,     "4",      8,
    # ANRHI missing - can NOT grade 0 or 2
    "Alkalosis",   7.5,    7.35,   NA,     NA,       9,
    "Alkalosis",   7.41,   7.35,   NA,     NA,       10,
    "Alkalosis",   7.4,    7.35,   NA,     NA,       11,
    "Alkalosis",   7.39,   7.35,   NA,     NA,       12,
    # AVAL missing cannot grade
    "Alkalosis",   NA,     1.1,    NA,     NA,       13,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_alkalo_daids, meta = atoxgr_criteria_daids)
})


### ALT, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 86: DAIDS ALT, High ----
test_that("derive_var_atoxgr_dir Test 86: DAIDS ALT, High", {
  local_atoxgr_test_data()
  exp_alti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",  30,     60,     NA,       1,
    NA_character_, 30,     60,     NA,       2,
    # ANRHI not missing
    "ALT, High",   601,    60,     "4",      3,
    "ALT, High",   600,    60,     "4",      4,
    "ALT, High",   599,    60,     "3",      5,
    "ALT, High",   300,    60,     "3",      6,
    "ALT, High",   299,    60,     "2",      7,
    "ALT, High",   150,    60,     "2",      8,
    "ALT, High",   149,    60,     "1",      9,
    "ALT, High",   76,     60,     "1",      10,
    "ALT, High",   75,     60,     "1",      11,
    "ALT, High",   74,     60,     "0",      12,
    # ANRHI missing cannot grade
    "ALT, High",   49,     NA,     NA,       13,
    # AVAL missing cannot grade
    "ALT, High",   NA,     60,     NA,       14,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_alti_daids, meta = atoxgr_criteria_daids)
})


### Amylase, High
### Grade 4: >= 5.0 x ULN
### Grade 3: 3.0 to < 5.0 x ULN
### Grade 2: 1.5 to < 3.0 x ULN
### Grade 1: 1.1 to < 1.5 x ULN

## Test 87: DAIDS Amylase, High ----
test_that("derive_var_atoxgr_dir Test 87: DAIDS Amylase, High", {
  local_atoxgr_test_data()
  exp_amyli_daids <- tibble::tribble(
    ~ATOXDSCH,         ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",      30,     60,     NA,       1,
    NA_character_,     30,     60,     NA,       2,
    # ANRHI not missing
    "Amylase, High",   301,    60,     "4",      3,
    "Amylase, High",   300,    60,     "4",      4,
    "Amylase, High",   299,    60,     "3",      5,
    "Amylase, High",   180,    60,     "3",      6,
    "Amylase, High",   179,    60,     "2",      7,
    "Amylase, High",   90,     60,     "2",      8,
    "Amylase, High",   89,     60,     "1",      9,
    "Amylase, High",   66,     60,     "1",      10,
    "Amylase, High",   65,     60,     "0",      11,
    # ANRHI missing cannot grade
    "Amylase, High",   65,     NA,     NA,       12,
    # AVAL missing cannot grade
    "Amylase, High",   NA,     60,     NA,       13,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_amyli_daids, meta = atoxgr_criteria_daids)
})

### AST, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 88: DAIDS AST, High ----
test_that("derive_var_atoxgr_dir Test 88: DAIDS AST, High", {
  local_atoxgr_test_data()
  exp_asti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",  30,     60,     NA,       1,
    NA_character_, 30,     60,     NA,       2,
    # ANRHI not missing
    "AST, High",   601,    60,     "4",      3,
    "AST, High",   600,    60,     "4",      4,
    "AST, High",   599,    60,     "3",      5,
    "AST, High",   300,    60,     "3",      6,
    "AST, High",   299,    60,     "2",      7,
    "AST, High",   150,    60,     "2",      8,
    "AST, High",   149,    60,     "1",      9,
    "AST, High",   76,     60,     "1",      10,
    "AST, High",   75,     60,     "1",      11,
    "AST, High",   74,     60,     "0",      12,
    # ANRHI missing cannot grade
    "AST, High",   49,     NA,     NA,       13,
    # AVAL missing cannot grade
    "AST, High",   NA,     60,     NA,       14,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_asti_daids, meta = atoxgr_criteria_daids)
})


### Bicarbonate, Low
### SI + cv unit is mmol/L
### Grade 4: < 8.0 mmol/L
### Grade 3: 8.0 -< 11.0 mmol/L
### Grade 2: 11.0 -< 16.0 mmol/L
### Grade 1: 16.0 mmol/L -< LLN


## Test 89a: DAIDS Bicarbonate, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 89a: DAIDS Bicarbonate, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_bicad_daids, meta = atoxgr_criteria_daids)
})

## Test 89b: DAIDS Bicarbonate, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 89b: DAIDS Bicarbonate, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_bicad_daids, meta = atoxgr_criteria_daids_uscv)
})


### Direct Bilirubin, High

### SI unit is umol/L
### CV unit is mg/dL
### 17.1 used as conversion from "mg/dL" to "umol/L"

### > 28 days of age

### Grade 4: > ULN

## SI unit

## CV unit

### <= 28 days of age

### Grade 4: > 2 mg/dL (> 34.2 umol/L)
### Grade 3: > 1.5 to <= 2 mg/dL (> 25.65 to <= 34.2 umol/L)
### Grade 2: > 1 to <= 1.5 mg/dL (> 17.1 to <= 25.65 umol/L)
### Grade 1: ULN to <= 1 mg/dL (ULN to <= 17.1 umol/L)

## SI unit

## CV unit

### add subjects with missing ADT or BRTHDT

## SI unit

## CV unit

### put all SI data together

### put all CV data together

## Test 90a: DAIDS Direct Bilirubin, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 90: DAIDS Direct Bilirubin, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_dbili_daids_si, meta = atoxgr_criteria_daids)
})

## Test 90b: DAIDS Direct Bilirubin, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 90b: DAIDS Direct Bilirubin, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_dbili_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Total Bilirubin, High

### > 28 days of age

### Grade 4: >= 5.0 x ULN
### Grade 3: 2.6 to < 5.0 x ULN
### Grade 2: 1.6 to < 2.6 x ULN
### Grade 1: 1.1 to < 1.6 x ULN


### make Age <= 28 all results NA for ATOXGRH

### make Age missing results NA for ATOXGRH



## Test 91: DAIDS Total Bilirubin, High ----
test_that("derive_var_atoxgr_dir Test 91: DAIDS Total Bilirubin, High", {
  local_atoxgr_test_data()
  test_high(expected = exp_tbili_daids, meta = atoxgr_criteria_daids)
})

### Calcium, High

## SI unit is mmol/L
### >= 7 days of age
### Grade 4: >= 3.38 mmol/L
### Grade 3: 3.13 -< 3.38 mmol/L
### Grade 2: 2.88 -< 3.13 mmol/L
### Grade 1: 2.65 -< 2.88 mmol/L

## CV unit is mg/dL
### >= 7 days of age
### Grade 4: >= 13.5 mg/dL
### Grade 3: 12.5 -< 13.5 mg/dL
### Grade 2: 11.5 -< 12.5 mg/dL
### Grade 1: 10.6 -< 11.5 mg/dL



### < 7 days of age
### SI unit is mmol/L
### Grade 4: >= 3.38 mmol/L
### Grade 3: 3.23 -< 3.38 mmol/L
### Grade 2: 3.1 -< 3.23 mmol/L
### Grade 1: 2.88 -< 3.1 mmol/L

## CV unit is mg/dL
### Grade 4: >= 13.5 mg/dL
### Grade 3: 12.9 -< 13.5 mg/dL
### Grade 2: 12.4 -< 12.9 mg/dL
### Grade 1: 11.5 -< 12.4 mg/dL



## create data with age missing for SI unit

## create data with age missing for CV unit

# put all SI data together

# put all CV data together
## Test 92a: DAIDS Calcium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 92a: DAIDS Calcium, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_daids_si, meta = atoxgr_criteria_daids)
})

## Test 92b: DAIDS Calcium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 92b: DAIDS Calcium, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calci_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Calcium (Ionized), High
### SI unit is mmol/L
### Grade 4: >= 1.8 mmol/L
### Grade 3: 1.6 -< 1.8 mmol/L
### Grade 2: 1.5 -< 1.6 mmol/L
### Grade 1: >ULN -< 1.5 mmol/L


### CV unit is mg/dL (4 * mmol/L)
### Grade 4: >= 7.2 mg/dL
### Grade 3: 6.4 -< 7.2 mg/dL
### Grade 2: 6.0 -< 6.4 mg/dL
### Grade 1: >ULN -< 6.0 mg/dL



## Test 93a: DAIDS Calcium (Ionized), High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 93a: DAIDS Calcium (Ionized), High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_daids_si, meta = atoxgr_criteria_daids)
})

## Test 93b: DAIDS Calcium (Ionized), High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 93b: DAIDS Calcium (Ionized), High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_calioni_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Calcium, Low

### >= 7 days of age

### SI unit is mmol/L
### Grade 4: < 1.53 mmol/L
### Grade 3: 1.53 -< 1.75 mmol/L
### Grade 2: 1.75 -< 1.95 mmol/L
### Grade 1: 1.95 -< 2.10 mmol/L

### CV unit is mg/dL
### Grade 4: < 6.1 mg/dL
### Grade 3: 6.1 -< 7.0 mg/dL
### Grade 2: 7.0 -< 7.8 mg/dL
### Grade 1: 7.8 -< 8.4 mg/dL



### < 7 days of age
### SI unit is mmol/L
### Grade 4: < 1.38 mmol/L
### Grade 3: 1.38 -< 1.5 mmol/L
### Grade 2: 1.5 -< 1.63 mmol/L
### Grade 1: 1.63 -< 1.88 mmol/L

### CV unit is mg/dL
### Grade 4: < 5.5 mg/dL
### Grade 3: 5.5 -< 6.0 mg/dL
### Grade 2: 6.0 -< 6.5 mg/dL
### Grade 1: 6.5 -< 7.5 mg/dL

### SI unit
### CV unit

## create SI records with age missing

## create CV records with age missing

## put all SI records together

## put all CV records together


## Test 94a: DAIDS Calcium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 94a: DAIDS Calcium, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_daids_si, meta = atoxgr_criteria_daids)
})

## Test 94b: DAIDS Calcium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 94b: DAIDS Calcium, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_calcd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Calcium (Ionized), Low
### SI unit is mmol/L
### Grade 4: <0.8 mmol/L
### Grade 3: 0.8 -< 0.9 mmol/L
### Grade 2: 0.9 -< 1.0 mmol/L
### Grade 1: 1.0 mmol/L -< LLN


### SI unit is mg/dL ( 4 * mmol/L)
### Grade 4: <3.2 mg/dL
### Grade 3: 3.2 -< 3.6 mg/dL
### Grade 2: 3.6 -< 4.0 mg/dL
### Grade 1: 4.0 mg/dL -< LLN



## Test 95a: DAIDS Calcium (Ionized), Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 95a: DAIDS Calcium (Ionized), Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_daids_si, meta = atoxgr_criteria_daids)
})

## Test 95b: DAIDS Calcium (Ionized), Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 95b: DAIDS Calcium (Ionized), Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_caliond_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Creatine Kinase, High
### Grade 4: >= 20.0 x ULN
### Grade 3: 10.0 -< 20.0 x ULN
### Grade 2: 6 -< 10.0 x ULN
### Grade 1: 3 -< 6 x ULN

## Test 96: DAIDS Creatine Kinase, High ----
test_that("derive_var_atoxgr_dir Test 96: DAIDS Creatine Kinase, High", {
  local_atoxgr_test_data()
  exp_cki_daids <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL, ~ANRHI, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",            10,    5,      NA_character_, NA,       1,
    NA_character_,           10,    5,      NA_character_, NA,       2,
    "Creatine Kinase, High", 100,   5,      NA_character_, "4",      3,
    "Creatine Kinase, High", 99,    5,      NA_character_, "3",      4,
    "Creatine Kinase, High", 50,    5,      NA_character_, "3",      5,
    "Creatine Kinase, High", 49,    5,      NA_character_, "2",      6,
    "Creatine Kinase, High", 30,    5,      NA_character_, "2",      7,
    "Creatine Kinase, High", 29,    5,      NA_character_, "1",      8,
    "Creatine Kinase, High", 15,    5,      NA_character_, "1",      9,
    "Creatine Kinase, High", 14,    5,      NA_character_, "0",      10,
    # ANRHI missing - cannot grade
    "Creatine Kinase, High", 4,     NA,     NA_character_, NA,       11,
    # AVAL missing cannot grade
    "Creatine Kinase, High", NA,    NA,     NA_character_, NA,       12,
  )

  test_high(expected = exp_cki_daids, meta = atoxgr_criteria_daids)
})


### Creatinine, High
### Grade 4: >= 3.5 x ULN or >= 2 X BASE
### Grade 3: >1.8 -< 3.5 x ULN or 1.5 -< 2 x BASE
### Grade 2: >1.3 - 1.8 x ULN or 1.3 - < 1.5 x BASE
### Grade 1: 1.1 - 1.3 x ULN

## Test 97: DAIDS Creatinine, High ----
test_that("derive_var_atoxgr_dir Test 97: DAIDS Creatinine, High", {
  local_atoxgr_test_data()
  exp_creati_daids <- tibble::tribble(
    ~ATOXDSCH,           ~AVAL, ~ANRHI, ~BASE, ~ATOXGRH, ~TESTNUM,
    "Not a term",        10,    10,     34,    NA,       1,
    NA_character_,       10,    10,     34,    NA,       2,
    "Creatinine, High",  35,    10,     34,    "4",      3,
    "Creatinine, High",  10,    10,     5,     "4",      4,
    "Creatinine, High",  34,    10,     34,    "3",      5,
    "Creatinine, High",  19,    10,     20,    "3",      6,
    "Creatinine, High",  9,     10,     5,     "3",      7,
    "Creatinine, High",  7.5,   10,     5,     "3",      8,
    "Creatinine, High",  18,    10,     34,    "2",      9,
    "Creatinine, High",  14,    10,     20,    "2",      10,
    "Creatinine, High",  7.4,   10,     5,     "2",      11,
    "Creatinine, High",  6.5,   10,     5,     "2",      12,
    "Creatinine, High",  13,    10,     34,    "1",      13,
    "Creatinine, High",  11,    10,     20,    "1",      14,
    "Creatinine, High",  10,    10,     20,    "0",      15,
    # ANRHI missing - cannot grade
    "Creatinine, High",  10,    NA,     20,    NA,       16,
    # AVAL missing cannot grade
    "Creatinine, High",  NA,    10,     20,    NA,       18,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_creati_daids, meta = atoxgr_criteria_daids)
})


### Glucose Fasting, High

### SI unit is mmol/L
### Grade 4: >= 27.75 mmol/L
### Grade 3: 13.89 -< 27.75 mmol/L
### Grade 2: 6.95 -< 13.89 mmol/L
### Grade 1: 6.11 -< 6.95 mmol/L

### CV unit is mg/dL
### Grade 4: >= 500 mg/dL
### Grade 3: >250 -< 500 mg/dL
### Grade 2: >125 - 250 mg/dL
### Grade 1: 110 - 125 mg/dL

## Test 98a:  DAIDS Glucose Fasting, High (SI unit)----
test_that("derive_var_atoxgr_dir Test 98a:  DAIDS Glucose Fasting, High (SI unit)", {
  local_atoxgr_test_data()
  exp_glucfi_daids_si <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
    "Not a term",            9.5,    "mg/L",   NA,       1,
    NA_character_,           9.5,    "mmol/L", NA,       2,
    "Glucose Fasting, High", 27.75,  "mmol/L", "4",      3,
    "Glucose Fasting, High", 27.74,  "mmol/L", "3",      4,
    "Glucose Fasting, High", 13.89,  "mmol/L", "3",      5,
    "Glucose Fasting, High", 13.88,  "mmol/L", "2",      6,
    "Glucose Fasting, High", 6.95,   "mmol/L", "2",      7,
    "Glucose Fasting, High", 6.94,   "mmol/L", "1",      8,
    "Glucose Fasting, High", 6.11,   "mmol/L", "1",      9,
    "Glucose Fasting, High", 6.1,    "mmol/L", "0",      10,
    # AVALU missing cannot grade
    "Glucose Fasting, High", 7,      NA,       NA,       11,
    # AVAL missing cannot grade
    "Glucose Fasting, High", NA,     "mmol/L", NA,       12,
  )

  test_high(expected = exp_glucfi_daids_si, meta = atoxgr_criteria_daids)
})

## Test 98b:  DAIDS Glucose Fasting, High (USCV unit)----
test_that("derive_var_atoxgr_dir Test 98b:  DAIDS Glucose Fasting, High (USCV unit)", {
  local_atoxgr_test_data()
  exp_glucfi_daids_cv <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
    "Not a term",            500,    "mg/L",   NA,       1,
    NA_character_,           500,    "mg/dL",  NA,       2,
    "Glucose Fasting, High", 500,    "mg/dL",  "4",      3,
    "Glucose Fasting, High", 499,    "mg/dL",  "3",      4,
    "Glucose Fasting, High", 251,    "mg/dL",  "3",      5,
    "Glucose Fasting, High", 250,    "mg/dL",  "2",      6,
    "Glucose Fasting, High", 126,    "mg/dL",  "2",      7,
    "Glucose Fasting, High", 125,    "mg/dL",  "1",      8,
    "Glucose Fasting, High", 110,    "mg/dL",  "1",      9,
    "Glucose Fasting, High", 109,    "mg/dL",  "0",      10,
    # AVALU missing cannot grade
    "Glucose Fasting, High", 7,      NA,       NA,       11,
    # AVAL missing cannot grade
    "Glucose Fasting, High", NA,     "mg/dL",  NA,       12,
  )

  test_high(expected = exp_glucfi_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Glucose Nonfasting, High

### SI unit is mmol/L
### Grade 4: >= 27.75 mmol/L
### Grade 3: 13.89 -< 27.75 mmol/L
### Grade 2: 8.89 -< 13.89 mmol/L
### Grade 1: 6.44 -< 8.89 mmol/L

### CV unit is mg/dL
### Grade 4: >= 500 mg/dL
### Grade 3: >250 -< 500 mg/dL
### Grade 2: >160 - 250 mg/dL
### Grade 1: 116 - 160 mg/dL

## Test 99a:  DAIDS Glucose Nonfasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 99a:  DAIDS Glucose Nonfasting, High (SI unit)", {
  local_atoxgr_test_data()
  exp_glucnfi_daids_si <- tibble::tribble(
    ~ATOXDSCH,                  ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
    "Not a term",               9.5,    "mg/L",   NA,       1,
    NA_character_,              9.5,    "mmol/L", NA,       2,
    "Glucose Nonfasting, High", 27.75,  "mmol/L", "4",      3,
    "Glucose Nonfasting, High", 27.74,  "mmol/L", "3",      4,
    "Glucose Nonfasting, High", 13.89,  "mmol/L", "3",      5,
    "Glucose Nonfasting, High", 13.88,  "mmol/L", "2",      6,
    "Glucose Nonfasting, High", 8.89,   "mmol/L", "2",      7,
    "Glucose Nonfasting, High", 8.88,   "mmol/L", "1",      8,
    "Glucose Nonfasting, High", 6.44,   "mmol/L", "1",      9,
    "Glucose Nonfasting, High", 6.43,   "mmol/L", "0",      10,
    # AVALU missing cannot grade
    "Glucose Nonfasting, High", 7,      NA,       NA,       11,
    # AVAL missing cannot grade
    "Glucose Nonfasting, High", NA,     "mmol/L", NA,       12,
  )

  test_high(expected = exp_glucnfi_daids_si, meta = atoxgr_criteria_daids)
})

## Test 99b:  DAIDS Glucose Nonfasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 99b:  DAIDS Glucose Nonfasting, High (USCV unit)", {
  local_atoxgr_test_data()
  exp_glucnfi_daids_cv <- tibble::tribble(
    ~ATOXDSCH,                   ~AVAL,  ~AVALU,  ~ATOXGRH,  ~TESTNUM,
    "Not a term",                500,    "mg/L",        NA,  1,
    NA_character_,               500,   "mg/dL",        NA,  2,
    "Glucose Nonfasting, High",  500,   "mg/dL",       "4",  3,
    "Glucose Nonfasting, High",  499,   "mg/dL",       "3",  4,
    "Glucose Nonfasting, High",  251,   "mg/dL",       "3",  5,
    "Glucose Nonfasting, High",  250,   "mg/dL",       "2",  6,
    "Glucose Nonfasting, High",  161,   "mg/dL",       "2",  7,
    "Glucose Nonfasting, High",  160,   "mg/dL",       "1",  8,
    "Glucose Nonfasting, High",  116,   "mg/dL",       "1",  9,
    "Glucose Nonfasting, High",  115,   "mg/dL",       "0",  10,
    # AVALU missing cannot grade
    "Glucose Nonfasting, High",  115,        NA,        NA,  11,
    # AVAL missing cannot grade
    "Glucose Nonfasting, High",  NA,    "mg/dL",        NA,  12,
  )

  test_high(expected = exp_glucnfi_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Glucose, Low

### >= 1 month of age

### SI unit is mmol/L
### Grade 4: < 1.67 mmol/L
### Grade 3: 1.67 -< 2.22 mmol/L
### Grade 2: 2.22 -< 3.05 mmol/L
### Grade 1: 3.05 -< 3.55 mmol/L

### CV unit is mg/dL
### Grade 4: < 30 mg/dL
### Grade 3: 30 -< 40 mg/dL
### Grade 2: 40 -< 55 mg/dL
### Grade 1: 55 - 64 mg/dL

### SI unit

### CV unit

### < 1 month of age

### SI unit is mmol/L
### Grade 4: < 1.67 mmol/L
### Grade 3: 1.67 -< 2.22 mmol/L
### Grade 2: 2.22 -< 2.78 mmol/L
### Grade 1: 2.78 -< 3.00 mmol/L

### CV unit is mg/dL
### Grade 4: < 30 mg/dL
### Grade 3: 30 -< 40 mg/dL
### Grade 2: 40 -< 50 mg/dL
### Grade 1: 50 - 54 mg/dL



## create SI records with age missing
## create CV records with age missing

## put all SI records together

## put all CV records together


## Test 100a:  DAIDS Glucose, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 100a:  DAIDS Glucose, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glucd_daids_si, meta = atoxgr_criteria_daids)
})

## Test 100b:  DAIDS Glucose, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 100b:  DAIDS Glucose, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_glucd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### Lactate, High

### Grade 2: >= 2.0 x ULN
### Grade 1: ULN -< 2.0 x ULN

## Test 101:  DAIDS Lactate, High ----
test_that("derive_var_atoxgr_dir Test 101:  DAIDS Lactate, High", {
  local_atoxgr_test_data()
  exp_lacti_daids <- tibble::tribble(
    ~ATOXDSCH,          ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",       105,    100,    NA,       1,
    NA_character_,      105,    100,    NA,       2,
    "Lactate, High",    200,    100,    "2",      3,
    "Lactate, High",    199,    100,    "1",      4,
    "Lactate, High",    100,    100,    "1",      5,
    "Lactate, High",    99,     100,    "0",      6,
    # ANRHI missing cannot grade
    "Lactate, High",    200,    NA,     NA,       7,
    # AVAL missing cannot grade
    "Lactate, High",    NA,     100,    NA,       8,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_lacti_daids, meta = atoxgr_criteria_daids)
})


### Lipase, High

### Grade 4: >= 5.0 x ULN
### Grade 3: 3.0 -< 5.0 x ULN
### Grade 2: 1.5 -< 3.0 x ULN
### Grade 1: 1.1 -< 1.5 x ULN

## Test 102:  DAIDS Lipase, High ----
test_that("derive_var_atoxgr_dir Test 102:  DAIDS Lipase, High", {
  local_atoxgr_test_data()
  exp_lipi_daids <- tibble::tribble(
    ~ATOXDSCH,          ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",       80,     100,    NA,       1,
    NA_character_,      60,     100,    NA,       2,
    "Lipase, High",     500,    100,    "4",      3,
    "Lipase, High",     499,    100,    "3",      4,
    "Lipase, High",     300,    100,    "3",      5,
    "Lipase, High",     299,    100,    "2",      6,
    "Lipase, High",     150,    100,    "2",      7,
    "Lipase, High",     149,    100,    "1",      8,
    "Lipase, High",     110,    100,    "1",      9,
    "Lipase, High",     109,    100,    "0",      10,
    # ANRHI missing cannot grade
    "Lipase, High",     200,    NA,     NA,       11,
    # AVAL missing cannot grade
    "Lipase, High",     NA,     100,    NA,       12,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_lipi_daids, meta = atoxgr_criteria_daids)
})


### Cholesterol, Fasting, High

### >= 18 years of age
### SI unit is mmol/L
### Grade 3: >= 7.77 mmol/L
### Grade 2: 6.19 -< 7.77 mmol/L
### Grade 1: 5.18 -< 6.19 mmol/L

### CV unit is mg/dL
### Grade 3: >= 300 mg/dL
### Grade 2: 240 -< 300 mg/dL
### Grade 1: 200 -< 240 mg/dL

### SI data

### CV data

### < 18 years of age
### SI unit is mmol/L
### Grade 3: >= 7.77 mmol/L
### Grade 2: 5.15 -< 7.77 mmol/L
### Grade 1: 4.4 -< 5.15 mmol/L

### CV unit is mg/dL
### Grade 3: >= 300 mg/dL
### Grade 2: 200 -< 300 mg/dL
### Grade 1: 170 -< 200 mg/dL



## create SI records with age missing

## create CV records with age missing

## put all SI data together

## put all CV data together

## Test 103a: DAIDS Cholesterol, Fasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 103a: DAIDS Cholesterol, Fasting, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_cholfi_daids_si, meta = atoxgr_criteria_daids)
})

## Test 103b: DAIDS Cholesterol, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 103b: DAIDS Cholesterol, Fasting, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_cholfi_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### LDL, Fasting, High

### >= 18 years of age

## SI unit is mmol/L
### Grade 3: >= 4.90 mmol/L
### Grade 2: 4.12 -< 4.90 mmol/L
### Grade 1: 3.17 -< 4.12 mmol/L

## CV unit is mg/dL
### Grade 3: >= 190 mg/dL
### Grade 2: 160 -< 190 mg/dL
### Grade 1: 130 -< 160 mg/dL



## SI unit is mmol/L
### > 2 to < 18 years of age
### Grade 3: >= 4.90 mmol/L
### Grade 2: 3.34 -< 4.90 mmol/L
### Grade 1: 2.85 -< 3.34 mmol/L

## CV unit is mg/dL
### Grade 3: >= 190 mg/dL
### Grade 2: 130 -< 190 mg/dL
### Grade 1: 110 -< 130 mg/dL

## SI data
## CV data

## create SI records with age mssing

## create CV records with age mssing

## create SI records with age <= 2 years

## create CV records with age <= 2 years

## put all SI records together

## put all CV records together

## Test 104a: DAIDS LDL, Fasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 104a: DAIDS LDL, Fasting, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_ldlfi_daids_si, meta = atoxgr_criteria_daids)
})

## Test 104b: DAIDS LDL, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 104b: DAIDS LDL, Fasting, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_ldlfi_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Triglycerides, Fasting, High

### SI unit is mmol/L
### Grade 4: > 11.4 mmol/L
### Grade 3: >5.7 - 11.4 mmol/L
### Grade 2: >3.42 - 5.7 mmol/L
### Grade 1: 1.71 - 3.42 mmol/L

### CV unit is mg/dL
### Grade 4: > 1000 mg/dL
### Grade 3: >500 - 1000 mg/dL
### Grade 2: >300 - 500 mg/dL
### Grade 1: 150 - 300 mg/dL

## Test 105a: DAIDS Triglycerides, Fasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 105a: DAIDS Triglycerides, Fasting, High (SI unit)", {
  local_atoxgr_test_data()
  exp_trigfi_daids_si <- tibble::tribble(
    ~ATOXDSCH,                       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
    "Not a term",                    1.5,    "mmol/L",  NA,       1,
    NA_character_,                   1.5,    "mmol/L",  NA,       2,
    "Triglycerides, Fasting, High",  11.5,   "mmol/L",  "4",      3,
    "Triglycerides, Fasting, High",  11.4,   "mmol/L",  "3",      4,
    "Triglycerides, Fasting, High",  5.8,    "mmol/L",  "3",      5,
    "Triglycerides, Fasting, High",  5.7,    "mmol/L",  "2",      6,
    "Triglycerides, Fasting, High",  3.43,   "mmol/L",  "2",      7,
    "Triglycerides, Fasting, High",  3.42,   "mmol/L",  "1",      8,
    "Triglycerides, Fasting, High",  1.71,   "mmol/L",  "1",      9,
    "Triglycerides, Fasting, High",  1.7,    "mmol/L",  "0",      10,
    # Unit missing cannot grade
    "Triglycerides, Fasting, High",  1.5,    NA,        NA,       11,
    # AVAL missing cannot grade
    "Triglycerides, Fasting, High",  NA,     "mmol/L",  NA,       12,
  )

  test_high(expected = exp_trigfi_daids_si, meta = atoxgr_criteria_daids)
})

## Test 105b: DAIDS Triglycerides, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 105b: DAIDS Triglycerides, Fasting, High (USCV unit)", {
  local_atoxgr_test_data()
  exp_trigfi_daids_cv <- tibble::tribble(
    ~ATOXDSCH,                       ~AVAL,   ~AVALU,    ~ATOXGRH, ~TESTNUM,
    "Not a term",                    1001,    "mg/dL",   NA,       1,
    NA_character_,                   1001,    "mg/dL",   NA,       2,
    "Triglycerides, Fasting, High",  1001,    "mg/dL",   "4",      3,
    "Triglycerides, Fasting, High",  1000,    "mg/dL",   "3",      4,
    "Triglycerides, Fasting, High",  501,     "mg/dL",   "3",      5,
    "Triglycerides, Fasting, High",  500,     "mg/dL",   "2",      6,
    "Triglycerides, Fasting, High",  301,     "mg/dL",   "2",      7,
    "Triglycerides, Fasting, High",  300,     "mg/dL",   "1",      8,
    "Triglycerides, Fasting, High",  150,     "mg/dL",   "1",      9,
    "Triglycerides, Fasting, High",  149,     "mg/dL",   "0",      10,
    # Unit missing cannot grade
    "Triglycerides, Fasting, High",  149,     NA,        NA,       11,
    # AVAL missing cannot grade
    "Triglycerides, Fasting, High",  NA,      "mg/dL",   NA,       12,
  )

  test_high(expected = exp_trigfi_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Magnesium, Low

### SI unit is mmol/L
### Grade 4: <0.3 mmol/L
### Grade 3: 0.3 -< 0.45 mmol/L
### Grade 2: 0.45 -< 0.6 mmol/L
### Grade 1: 0.6 -< 0.7 mmol/L

### CV unit is mg/dL
### Conversion factor 2.431 used to convert mmol/L to mg/dL


## CV unit using conversions factor 2.431

## Test 106a: DAIDS Magnesium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 106a: DAIDS Magnesium, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magd_daids_si, meta = atoxgr_criteria_daids)
})

## Test 106b: DAIDS Magnesium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 106b: DAIDS Magnesium, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_magd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Phosphate, Low

### > 14 years of age

### SI unit is mmol/L
### Grade 4: < 0.32 mmol/L
### Grade 3: 0.32 to < 0.45 mmol/L
### Grade 2: 0.45 to < 0.65 mmol/L
### Grade 1: 0.65 mmol/L to < LLN mmol/L

### CV unit is mg/dL
### Grade 4: < 1 mg/dL
### Grade 3: 1 to < 1.4 mg/dL
### Grade 2: 1.4 to < 2 mg/dL
### Grade 1: 2 to < LLN mg/dL



### 1 to 14 years of age

### SI unit is mmol/L
### Grade 4: < 0.48 mmol/L
### Grade 3: 0.48 to < 0.81 mmol/L
### Grade 2: 0.81 to < 0.97 mmol/L
### Grade 1: 0.97 to <1.13 mmol/L

### CV unit is mg/dL
### Grade 4: < 1.5 mg/dL
### Grade 3: 1.5 to < 2.5 mg/dL
### Grade 2: 2.5 to < 3 mg/dL
### Grade 1: 3 to < 3.5 mg/dL



### < 1 year of age

### SI unit is mmol/L
### Grade 4: < 0.48 mmol/L
### Grade 3: 0.48 to < 0.81 mmol/L
### Grade 2: 0.81 to < 1.13 mmol/L
### Grade 1: 1.13 to < 1.45 mmol/L

### CV unit is mg/dL
### Grade 4: < 1.5 mg/dL
### Grade 3: 1.5 to < 2.5 mg/dL
### Grade 2: 2.5 to < 3.5 mg/dL
### Grade 1: 3.5 to < 4.5 mg/dL




# Set lab date or birth date to missing for SI records
# Set lab date or birth date to missing for CV records

## put all SI records together

## put all CV records together

## Test 107a: DAIDS Phosphate, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 107a: DAIDS Phosphate, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_phosd_daids_si, meta = atoxgr_criteria_daids)
})

## Test 107b: DAIDS Phosphate, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 107b: DAIDS Phosphate, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_phosd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Potassium, High

### SI unit + CV is mmol/L
### Grade 4: >= 7 mmol/L
### Grade 3: 6.5 -< 7 mmol/L
### Grade 2: 6 -< 6.5 mmol/L
### Grade 1: 5.6 -< 6 mmol/L


## Test 108a: DAIDS Potassium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 108a: DAIDS Potassium, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_poti_daids, meta = atoxgr_criteria_daids)
})

## Test 108b: DAIDS Potassium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 108b: DAIDS Potassium, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_poti_daids, meta = atoxgr_criteria_daids_uscv)
})

### Potassium, Low

### SI unit + CV is mmol/L
### Grade 4: <2 mmol/L
### Grade 3: 2 -< 2.5 mmol/L
### Grade 2: 2.5 -< 3 mmol/L
### Grade 1: 3 -< 3.4 mmol/L


## Test 109a: DAIDS Potassium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 109a: DAIDS Potassium, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_potd_daids, meta = atoxgr_criteria_daids)
})

## Test 109b: DAIDS Potassium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 109b: DAIDS Potassium, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_potd_daids, meta = atoxgr_criteria_daids_uscv)
})

### Sodium, High

### SI unit + CV is mmol/L
### Grade 4: >= 160 mmol/L
### Grade 3: 154 -< 160 mmol/L
### Grade 2: 150 -< 154 mmol/L
### Grade 1: 146 -< 150 mmol/L


## Test 110a: DAIDS Sodium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 110a: DAIDS Sodium, High (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_sodi_daids, meta = atoxgr_criteria_daids)
})

## Test 110b: DAIDS Sodium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 110b: DAIDS Sodium, High (USCV unit)", {
  local_atoxgr_test_data()
  test_high(exp_sodi_daids, meta = atoxgr_criteria_daids_uscv)
})

### Sodium, Low

### SI unit + CV is mmol/L
### Grade 4: <= 120 mmol/L
### Grade 3: >120 -< 125 mmol/L
### Grade 2: 125 -< 130 mmol/L
### Grade 1: 130 -< 135 mmol/L


## Test 111a: DAIDS Sodium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 111a: DAIDS Sodium, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_sodd_daids, meta = atoxgr_criteria_daids)
})

## Test 111b: DAIDS Sodium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 111b: DAIDS Sodium, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_sodd_daids, meta = atoxgr_criteria_daids_uscv)
})
### Uric Acid, High

### SI unit is umol/L
### Grade 4: >= 890 umol/L (0.89 mmol/L)
### Grade 3: 710 -< 890 umol/L (0.71 - <0.89 mmol/L)
### Grade 2: 590 -< 710 umol/L (0.59 - <0.71 mmol/L)
### Grade 1: 450 -< 590 umol/L (0.45 - <0.59 mmol/L)

### CV unit is mg/dL
### Grade 4: >= 15 mg/dL
### Grade 3: 12 -< 15 mg/dL
### Grade 2: 10 -< 12 mg/dL
### Grade 1: 7.5 -< 10 mg/dL

## Test 112a: DAIDS Uric Acid, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 112a: DAIDS Uric Acid, High (SI unit)", {
  local_atoxgr_test_data()
  exp_urici_daids <- tibble::tribble(
    ~ATOXDSCH,          ~AVAL, ~AVALU,    ~ATOXGRH, ~TESTNUM,
    "Not a term",       591,   "umol/L",  NA,       1,
    NA_character_,      591,   "umol/L",  NA,       2,
    "Uric Acid, High",  890,   "umol/L",  "4",      3,
    "Uric Acid, High",  889,   "umol/L",  "3",      4,
    "Uric Acid, High",  710,   "umol/L",  "3",      5,
    "Uric Acid, High",  709,   "umol/L",  "2",      6,
    "Uric Acid, High",  590,   "umol/L",  "2",      7,
    "Uric Acid, High",  589,   "umol/L",  "1",      8,
    "Uric Acid, High",  450,   "umol/L",  "1",      9,
    "Uric Acid, High",  449,   "umol/L",  "0",      10,
    # Unit missing cannot grade
    "Uric Acid, High",  200,   NA,        NA,       11,
    # AVAL missing cannot grade
    "Uric Acid, High",  NA,    "umol/L",  NA,       12,
  )

  test_high(expected = exp_urici_daids, meta = atoxgr_criteria_daids)
})

## Test 112b: DAIDS Uric Acid, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 112b: DAIDS Uric Acid, High (USCV unit)", {
  local_atoxgr_test_data()
  exp_uric_daids <- tibble::tribble(
    ~ATOXDSCH,         ~AVAL,  ~AVALU,  ~ATOXGRH, ~TESTNUM,
    "Not a term",      15,    "mg/dL",        NA, 1,
    NA_character_,     15,    "mg/dL",        NA, 2,
    "Uric Acid, High", 15,    "mg/dL",       "4", 3,
    "Uric Acid, High", 14.9,  "mg/dL",       "3", 4,
    "Uric Acid, High", 12,    "mg/dL",       "3", 5,
    "Uric Acid, High", 11.9,  "mg/dL",       "2", 6,
    "Uric Acid, High", 10,    "mg/dL",       "2", 7,
    "Uric Acid, High", 9.9,   "mg/dL",       "1", 8,
    "Uric Acid, High", 7.5,   "mg/dL",       "1", 9,
    "Uric Acid, High", 7.4,   "mg/dL",       "0", 10,
    # Unit missing cannot grade
    "Uric Acid, High", 200,        NA,        NA, 11,
    # AVAL missing cannot grade
    "Uric Acid, High", NA,    "mg/dL",        NA, 12,
  )

  test_high(expected = exp_uric_daids, meta = atoxgr_criteria_daids_uscv)
})

### Absolute CD4+ Count, Low

### > 5 years of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 0.1 10^9/L
### Grade 3: 0.1 to < 0.2 10^9/L
### Grade 2: 0.2 to < 0.3 10^9/L
### Grade 1: 0.3 to < 0.4 10^9/L


# no criteria for age <= 5 years set grade to missing

# add missing ADT and BRTHDT


## Test 113a: DAIDS Absolute CD4 Count, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 113a: DAIDS Absolute CD4 Count, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4d_daids_si, meta = atoxgr_criteria_daids)
})


## Test 113b: DAIDS Absolute CD4 Count, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 113b: DAIDS Absolute CD4 Count, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_cd4d_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Absolute Lymphocyte Count, Low

### > 5 years of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 0.35 10^9/L
### Grade 3: 0.35 to < 0.5 10^9/L
### Grade 2: 0.5 to < 0.6 10^9/L
### Grade 1: 0.6 to < 0.65 10^9/L


# no criteria for age <= 5 years set grade to missing

# add missing ADT and BRTHDT


## Test 114a: DAIDS Absolute Lymphocyte Count, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 114a: DAIDS Absolute Lymphocyte Count, Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymphd_daids_si, meta = atoxgr_criteria_daids)
})


## Test 114b: DAIDS Absolute Lymphocyte Count, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 114b: DAIDS Absolute Lymphocyte Count, Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymphd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


## Test 114c: DAIDS Absolute Lymph. Count, Low (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 114c: DAIDS Absolute Lymph. Count, Low (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_lymphd_daids_cv2, meta = atoxgr_criteria_daids_uscv)
})

### Absolute Neutrophil Count (ANC), Low

### > 7 days of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 0.4 10^9/L
### Grade 3: 0.4 to < 0.6 10^9/L
### Grade 2: 0.6 to < 0.8 10^9/L
### Grade 1: 0.8 to 1 10^9/L


### 2 to 7 days of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 0.75 10^9/L
### Grade 3: 0.75 to < 1.0 10^9/L
### Grade 2: 1.0 to < 1.25 10^9/L
### Grade 1: 1.25 to 1.5 10^9/L


### <= 1 day of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 1.50 10^9/L
### Grade 3: 1.50 to < 3.0 10^9/L
### Grade 2: 3.0 to < 4.0 10^9/L
### Grade 1: 4.0 to 5.0 10^9/L



# Set lab date/birth date to missing


## Test 115a: DAIDS ANC Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 115a: DAIDS ANC Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_ancd_daids_si, meta = atoxgr_criteria_daids)
})


## Test 115b: DAIDS ANC Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 115b: DAIDS ANC Low (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_ancd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


## Test 115c: DAIDS ANC Low (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 115c: DAIDS ANC Low (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_ancd_daids_cv2, meta = atoxgr_criteria_daids_uscv)
})


### Fibrinogen Decreased
### SI unit is g/L
### Grade 4: <0.5 g/L OR < 0.25 x LLN
### Grade 3: 0.5 to <0.75 g/L OR 0.25 to < 0.50 x LLN
### Grade 2: 0.75 to <1 g/L OR >= 0.50 to < 0.75 x LLN
### Grade 1: 1 to < 2 g/L OR 0.75 to < 1.00 x LLN

### CV unit is mg/dL (= 100 * g/L)


## Test 116a: DAIDS Fibrinogen Decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 116a: DAIDS Fibrinogen Decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_fibd_daids_si, meta = atoxgr_criteria_daids)
})

## CV unit is mg/dL which is 100 * SI unit of g/L

## Test 116b: DAIDS Fibrinogen Decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 116b: DAIDS Fibrinogen Decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_fibd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})

### Hemoglobin, Low
### SI unit is g/L
### CV unit g/dL ( = g/L / 10)

### >= 13 years of age (male only)

### Grade 4: < 70 g/L
### Grade 3: 70 to < 90 g/L
### Grade 2: 90 to < 100 g/L
### Grade 1: 100 to 109 g/L


### >= 13 years of age (female only)

### Grade 4: < 65 g/L
### Grade 3: 65 to < 85 g/L
### Grade 2: 85 to < 95 g/L
### Grade 1: 95 to 104 g/L



### 57 days to < 13 years of age (male and female)

### Grade 4: < 65 g/L
### Grade 3: 65 to < 85 g/L
### Grade 2: 85 to < 95 g/L
### Grade 1: 95 to 104 g/L


### 36 to <= 56 days of age (male and female)

### Grade 4: < 60 g/L
### Grade 3: 60 to < 70 g/L
### Grade 2: 70 to < 85 g/L
### Grade 1: 85 to 96 g/L



### 22 to <= 35 days of age (male and female)

### Grade 4: < 67 g/L
### Grade 3: 67 to < 80 g/L
### Grade 2: 80 to < 95 g/L
### Grade 1: 95 to 110 g/L



### 8 to <= 21 days of age (male and female)

### Grade 4: < 80 g/L
### Grade 3: 80 to < 90 g/L
### Grade 2: 90 to < 110 g/L
### Grade 1: 110 to 130 g/L


### <= 7 days of age (male and female)

### Grade 4: < 90 g/L
### Grade 3: 90 to < 100 g/L
### Grade 2: 100 to < 130 g/L
### Grade 1: 130 to 140 g/L



# Set lab date to missing fo each type, ie SEX is M, F or missing


## Test 117a: DAIDS HGB Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 117a: DAIDS HGB Low (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_hgbd_daids, meta = atoxgr_criteria_daids)
})

## Test 117b: DAIDS HGB Low (CV unit) ----
test_that("derive_var_atoxgr_dir Test 117b: DAIDS HGB Low (CV unit)", {
  local_atoxgr_test_data()
  exp_hgbd_daids_cv <- exp_hgbd_daids %>%
    mutate(
      AVAL = AVAL / 10,
      AVALU = if_else(!is.na(AVALU) & str_to_upper(AVALU) == "G/L", "g/dL", NA_character_)
    )

  test_low(expected = exp_hgbd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


### INR, high

### Grade 4: >=3 x ULN
### Grade 3: 2 to <3 x ULN
### Grade 2: 1.5 to < 2 x ULN
### Grade 1: 1.1 to < 1.5 x ULN

## Test 118: DAIDS INR High ----
test_that("derive_var_atoxgr_dir Test 118: DAIDS INR High", {
  local_atoxgr_test_data()
  exp_inri_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",  80,     80,     NA,       1,
    NA_character_, 60,     80,     NA,       2,
    "INR, High",   240,    80,     "4",      3,
    "INR, High",   239,    80,     "3",      4,
    "INR, High",   160,    80,     "3",      5,
    "INR, High",   159,    80,     "2",      6,
    "INR, High",   120,    80,     "2",      7,
    "INR, High",   119,    80,     "1",      8,
    "INR, High",   88,     80,     "1",      9,
    "INR, High",   87,     80,     "0",      10,
    # ANRHI missing - cannot grade
    "INR, High",   100,    NA,     NA,       11,
    # AVAL missing cannot grade
    "INR, High",   NA,     80,     NA,       12,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_inri_daids, meta = atoxgr_criteria_daids)
})


### Methemoglobin
### SI + CV unit is %

### Grade 4: >=20.0%
### Grade 3: 15 to < 20%
### Grade 2: 10 to < 15%
### Grade 1: 5 to <10%


## Test 119a: DAIDS Methemoglobin (SI unit) ----
test_that("derive_var_atoxgr_dir Test 119a: DAIDS Methemoglobin (SI unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_methi_daids, meta = atoxgr_criteria_daids)
})

## Test 119b: DAIDS Methemoglobin (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 119b: DAIDS Methemoglobin (USCV unit)", {
  local_atoxgr_test_data()
  test_high(expected = exp_methi_daids, meta = atoxgr_criteria_daids_uscv)
})

### PTT, high

### Grade 4: >=3 x ULN
### Grade 3: 2.33 to <3 x ULN
### Grade 2: 1.66 to < 2.33 x ULN
### Grade 1: 1.1 to < 1.66 x ULN

## Test 120: DAIDS PTT High ----
test_that("derive_var_atoxgr_dir Test 120: DAIDS PTT High", {
  local_atoxgr_test_data()
  exp_ptti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL, ~ANRHI, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",  80,    80,     NA_character_, NA,       1,
    NA_character_, 60,    80,     NA_character_, NA,       2,
    "PTT, High",   240,   80,     NA_character_, "4",      3,
    "PTT, High",   239,   80,     NA_character_, "3",      4,
    "PTT, High",   186.4, 80,     NA_character_, "3",      5,
    "PTT, High",   186.3, 80,     NA_character_, "2",      6,
    "PTT, High",   132.8, 80,     NA_character_, "2",      7,
    "PTT, High",   132.7, 80,     NA_character_, "1",      8,
    "PTT, High",   88,    80,     NA_character_, "1",      9,
    "PTT, High",   87,    80,     NA_character_, "0",      10,
    # ANRHI missing - cannot grade
    "PTT, High",   100,   NA,     NA_character_, NA,       11,
    # AVAL missing cannot grade
    "PTT, High",   NA,    80,     NA_character_, NA,       12,
  )

  test_high(expected = exp_ptti_daids, meta = atoxgr_criteria_daids)
})


### Platelets, Decreased
### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <25 x 10e9 /L
### Grade 3: 25 to <50  x 10e9 /L
### Grade 2: 50 to <100 - x 10e9
### Grade 1: 100 - 125 x 10e9 /L


## Test 121a: DAIDS Platelets decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 121a: DAIDS Platelets decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plated_daids_si, meta = atoxgr_criteria_daids)
})


## Test 121b: DAIDS Platelets decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 121b: DAIDS Platelets decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plated_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


## Test 121c: DAIDS Platelets decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 121c: DAIDS Platelets decreased (legacy USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_plated_daids_cv2, meta = atoxgr_criteria_daids_uscv)
})

### PT, high

### Grade 4: >=3 x ULN
### Grade 3: 1.5 - <3 x ULN
### Grade 2: 1.25 - <1.5 x ULN
### Grade 1: 1.1 - <1.25 x ULN

## Test 122: DAIDS PT High ----
test_that("derive_var_atoxgr_dir Test 122: DAIDS PT High", {
  local_atoxgr_test_data()
  exp_pti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI,  ~AVALU,         ~ATOXGRH, ~TESTNUM,
    "Not a term",  80,     100,     NA_character_,  NA,       1,
    NA_character_, 60,     100,     NA_character_,  NA,       2,
    "PT, High",    300,    100,     NA_character_,  "4",      3,
    "PT, High",    299,    100,     NA_character_,  "3",      4,
    "PT, High",    150,    100,     NA_character_,  "3",      5,
    "PT, High",    149,    100,     NA_character_,  "2",      6,
    "PT, High",    125,    100,     NA_character_,  "2",      7,
    "PT, High",    124,    100,     NA_character_,  "1",      8,
    "PT, High",    110,    100,     NA_character_,  "1",      9,
    "PT, High",    109,    100,     NA_character_,  "0",      10,
    # ANRHI missing - cannot grade
    "PT, High",    100,    NA,      NA_character_,  NA,       11,
    # AVAL missing cannot grade
    "PT, High",    NA,     100,     NA_character_,  NA,       12,
  )

  test_high(expected = exp_pti_daids, meta = atoxgr_criteria_daids)
})

### White blood cell decreased (> 7 days of age)

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <1 x 10e9/L
### Grade 3: 1 to 1.499 x 10e9/L
### Grade 2: 1.5 to 1.999 x 10e9/L
### Grade 1: 2 to 2.499 x 10e9/L


### White blood cell decreased (<= 7 days of age)

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <2.500 x 10e9/L
### Grade 3: 2.5 to 3.999 x 10e9/L
### Grade 2:   4 to 5.499 x 10e9/L
### Grade 1: 5.5 to 6.999 x 10e9/L




## Test 123a: DAIDS White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 123a: DAIDS White blood cell decreased (SI unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_daids_si, meta = atoxgr_criteria_daids)
})


## Test 123b: DAIDS White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 123b: DAIDS White blood cell decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_daids_cv, meta = atoxgr_criteria_daids_uscv)
})


## Test 123c: DAIDS White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 123c: DAIDS White blood cell decreased (USCV unit)", {
  local_atoxgr_test_data()
  test_low(expected = exp_wbcd_daids_cv2, meta = atoxgr_criteria_daids_uscv)
})

# NCICTCAEv6 new grading----

### Blood lactate dehydrogenase increased

### Grade 1: >ULN

## Test 124: CTCAEv6 LDH (SI unit) ----
test_that("derive_var_atoxgr_dir Test 124: CTCAEv6 LDH increased (SI unit)", {
  local_atoxgr_test_data()
  exp_ldhi_si <- tibble::tribble(
    ~ATOXDSCH,                               ~AVAL, ~ANRHI,      ~ATOXGRH,  ~TESTNUM,
    "Not a term",                              1.1,      1, NA_character_,         1,
    NA_character_,                             1.1,      1, NA_character_,         2,
    "Blood lactate dehydrogenase increased",   1.1,      1,           "1",         3,
    "Blood lactate dehydrogenase increased",     1,      1,           "0",         4,
    "Blood lactate dehydrogenase increased",    NA,      1, NA_character_,         5,
    "Blood lactate dehydrogenase increased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_ldhi_si, meta = atoxgr_criteria_ctcv6)
})

## Test 125: CTCAEv6 HDL decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 125: CTCAEv6 HDL decreased (SI unit)", {
  local_atoxgr_test_data()
  exp_hdld_si <- tibble::tribble(
    ~ATOXDSCL,       ~AVAL, ~ANRLO,      ~ATOXGRL,  ~TESTNUM,
    "Not a term",      0.9,      1, NA_character_,         1,
    NA_character_,     0.9,      1, NA_character_,         2,
    "HDL decreased",   0.9,      1,           "1",         3,
    "HDL decreased",     1,      1,           "0",         4,
    "HDL decreased",    NA,      1, NA_character_,         5,
    "HDL decreased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  test_low(expected = exp_hdld_si, meta = atoxgr_criteria_ctcv6)
})

## Test 126: CTCAEv6 LDL increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 126: CTCAEv6 LDL increased (SI unit)", {
  local_atoxgr_test_data()
  exp_ldli_si <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRHI,      ~ATOXGRH,  ~TESTNUM,
    "Not a term",      1.1,      1, NA_character_,         1,
    NA_character_,     1.1,      1, NA_character_,         2,
    "LDL increased",   1.1,      1,           "1",         3,
    "LDL increased",     1,      1,           "0",         4,
    "LDL increased",    NA,      1, NA_character_,         5,
    "LDL increased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_ldli_si, meta = atoxgr_criteria_ctcv6)
})

### Acidosis
### GRADE_3: pH is <7.3
### GRADE_1: pH is <normal, but >=7.3

## Test 127: CTCAEv6 Acidosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 127: CTCAEv6 Acidosis (SI unit)", {
  local_atoxgr_test_data()
  exp_acido <- tibble::tribble(
    ~ATOXDSCL,     ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRL, ~TESTNUM,
    "Not a term",  7.3,    7.35,   7.4,    NA,       1,
    NA_character_, 7.3,    7.35,   7.4,    NA,       2,
    # ANRLO not missing
    "Acidosis",    7.29,   7.35,   7.4,    "3",      3,
    "Acidosis",    7.3,    7.35,   7.4,    "1",      4,
    "Acidosis",    7.34,   7.35,   7.4,    "1",      5,
    "Acidosis",    7.35,   7.35,   7.4,    "0",      6,
    "Acidosis",    7.36,   7.35,   7.4,    "0",      7,
    # ANRLO missing - can grade 4
    "Acidosis",    7.29,   NA,     7.4,    "3",      8,
    # ANRLO missing - can NOT grade 0 or 1
    "Acidosis",    7.3,    NA,     7.4,    NA,       9,
    "Acidosis",    7.34,   NA,     7.4,    NA,       10,
    "Acidosis",    7.35,   NA,     7.4,    NA,       11,
    "Acidosis",    7.36,   NA,     7.4,    NA,       12,
    # AVAL missing cannot grade
    "Acidosis",    NA,     1.1,    1.4,    NA,       13,
  ) %>%
    mutate(AVALU = NA_character_)

  test_low(expected = exp_acido, meta = atoxgr_criteria_ctcv6)
})


### Acidosis
### GRADE_3: pH is >7.5
### GRADE_1: pH is >normal, but <=7.5

## Test 128: CTCAEv6 Alkalosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 128: CTCAEv6 Alkalosis (SI unit)", {
  local_atoxgr_test_data()
  exp_alkalo <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
    "Not a term",  7.3,    7.35,   7.4,    NA,       1,
    NA_character_, 7.3,    7.35,   7.4,    NA,       2,
    # ANRHI not missing
    "Alkalosis",   7.51,   7.35,   7.4,    "3",      3,
    "Alkalosis",   7.5,    7.35,   7.4,    "1",      4,
    "Alkalosis",   7.41,   7.35,   7.4,    "1",      5,
    "Alkalosis",   7.4,    7.35,   7.4,    "0",      6,
    "Alkalosis",   7.39,   7.35,   7.4,    "0",      7,
    # ANRHI missing - can grade 3
    "Alkalosis",   7.51,   7.35,   NA,     "3",      8,
    # ANRHI missing - can NOT grade 0 or 1
    "Alkalosis",   7.5,    7.35,   NA,     NA,       9,
    "Alkalosis",   7.41,   7.35,   NA,     NA,       10,
    "Alkalosis",   7.4,    7.35,   NA,     NA,       11,
    "Alkalosis",   7.39,   7.35,   NA,     NA,       12,
    # AVAL missing cannot grade
    "Alkalosis",   NA,     1.1,    NA,     NA,       13,
  ) %>%
    mutate(AVALU = NA_character_)

  test_high(expected = exp_alkalo, meta = atoxgr_criteria_ctcv6)
})

### Creatinine clearance decreased
### GRADE_3: 10 - 25 ml/min
### GRADE_2: 26 - 49 ml/min

## Test 129: CTCAEv6 Creat. clear. dec. (SI unit) ----
test_that("derive_var_atoxgr_dir Test 129: CTCAEv6 Creat. clear. dec. (SI unit)", {
  local_atoxgr_test_data()
  exp_crcl_d <- tibble::tribble(
    ~ATOXDSCL,                         ~AVAL,  ~ATOXGRL,  ~AVALU,    ~TESTNUM,
    "Not a term",                      27,     NA,        "mL/min",  1,
    NA_character_,                     27,     NA,        "mL/min",  2,
    # AVAL < 10 - should be missing grade
    "Creatinine clearance decreased",  9,      NA,        "mL/min",  3,
    # get grade 3 and 2
    "Creatinine clearance decreased",  10,     "3",       "mL/min",  4,
    "Creatinine clearance decreased",  25,     "3",       "mL/min",  5,
    "Creatinine clearance decreased",  26,     "2",       "mL/min",  6,
    "Creatinine clearance decreased",  49,     "2",       "mL/min",  7,
    "Creatinine clearance decreased",  50,     "0",       "mL/min",  8,
    # AVAL missing cannot grade
    "Creatinine clearance decreased",  NA,     NA,        "mL/min",  13,
    # unit incorrect cannot grade
    "Creatinine clearance decreased",  50,     NA,        "L/min",   14,
  )

  test_low(expected = exp_crcl_d, meta = atoxgr_criteria_ctcv6)
})

## Test when deprecated abnormal_indicator used - should map to high_indicator
## Test 130: CTCAEv6  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 130: CTCAEv6  Blood bilirubin increased", {
  local_atoxgr_test_data()
  input_bili_ctcv6 <- exp_bili_ctcv6 %>%
    select(-ATOXGRH)

  expect_snapshot(
    actual_bili_ctcv6 <- derive_var_atoxgr_dir(
      input_bili_ctcv6,
      new_var = ATOXGRH,
      meta_criteria = atoxgr_criteria_ctcv6,
      tox_description_var = ATOXDSCH,
      criteria_direction = "H",
      abnormal_indicator = "HIGH",
      get_unit_expr = AVALU
    )
  )

  expect_dfs_equal(
    base = exp_bili_ctcv6,
    compare = actual_bili_ctcv6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
})

## Test when high_indicator not defined - should map to high_indicator
## Test 131: CTCAEv6  Blood bili incr. high_indicator not defined ----
test_that("derive_var_atoxgr_dir Test 131: CTCAEv6  Blood bili incr. high_indicator not defined", {
  local_atoxgr_test_data()
  input_bili_ctcv6 <- exp_bili_ctcv6 %>%
    select(-ATOXGRH)

  expect_error(
    actual_bili_ctcv6 <- derive_var_atoxgr_dir(
      input_bili_ctcv6,
      new_var = ATOXGRH,
      meta_criteria = atoxgr_criteria_ctcv6,
      tox_description_var = ATOXDSCH,
      criteria_direction = "H",
      high_indicator = NULL,
      get_unit_expr = AVALU
    ),
    class = "assert_character_vector"
  )
})

## Test 132: CTCAEv6 Creatinine incr. low_indicator not defined ----
test_that("derive_var_atoxgr_dir Test 132: CTCAEv6 Creatinine incr. low_indicator not defined", {
  local_atoxgr_test_data()
  exp_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~BNRIND,        ~ATOXGRH,
    "Creatinine increased",  241,    40,    40,     NA_character_,  "NORMAL",       "4",
  )

  input_creatn <- exp_creatn %>%
    select(-ATOXGRH)

  expect_error(
    actual_creatn <- derive_var_atoxgr_dir(
      input_creatn,
      new_var = ATOXGRH,
      meta_criteria = atoxgr_criteria_ctcv6,
      tox_description_var = ATOXDSCH,
      criteria_direction = "H",
      high_indicator = "HIGH",
      low_indicator = NULL,
      get_unit_expr = AVALU
    ),
    class = "assert_character_vector"
  )
})

## Test 133: ERROR when CRITERIA different with same TERM and UNIT ----
test_that("derive_var_atoxgr_dir Test 133: ERROR when CRITERIA different with same TERM and UNIT", {
  local_atoxgr_test_data()
  atoxgr_criteria_ctcv6 <- bind_rows(
    admiral::atoxgr_criteria_ctcv6,
    admiral::atoxgr_criteria_ctcv6 %>%
      filter(TERM == "Creatinine increased") %>%
      mutate(GRADE_CRITERIA_CODE = paste0(GRADE_CRITERIA_CODE, ","))
  )

  exp_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~BNRIND,        ~ATOXGRH,
    "Creatinine increased",  241,    40,    40,     NA_character_,  "NORMAL",       "4",
  )

  input_creatn <- exp_creatn %>%
    select(-ATOXGRH)

  expect_error(
    actual_creatn <- derive_var_atoxgr_dir(
      input_creatn,
      new_var = ATOXGRH,
      meta_criteria = atoxgr_criteria_ctcv6,
      tox_description_var = ATOXDSCH,
      criteria_direction = "H",
      high_indicator = "HIGH",
      low_indicator = "LOW",
      get_unit_expr = AVALU
    ),
    NULL
  )
})
