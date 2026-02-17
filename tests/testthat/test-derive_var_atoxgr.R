## Test 1: ATOXGR cannot be graded ----
test_that("derive_var_atoxgr Test 1: ATOXGR cannot be graded", {
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

## Blood and lymphatic system disorders

### Anemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: <80 g/L
### Grade 2: <100 - 80g/L
### Grade 1: <LLN - 100 g/L

expected_anemia_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",   80,     120,    200,    "G/L",  NA,       1,
  NA_character_,  60,     50,     100,    "G/L",  NA,       2,
  "ANEMIA",       79,     140,    NA,     "G/L",  "3",      3,
  "ANEMIA",       80,     140,    NA,     "G/L",  "2",      4,
  "Anemia",       99,     140,    NA,     "G/L",  "2",      5,
  "Anemia",       100,    140,    NA,     "G/L",  "1",      6,
  # wrong UNIT - GRADE should be missing
  "anemia",       100,    140,    NA,     "G/cL", NA,       7,
  "Anemia",       139,    140,    NA,     "G/L",  "1",      8,
  "ANEMIA",       140,    140,    NA,     "g/L",  "0",      9,
  # ANRLO missing AVAL not grade 2 or 3 - cannot grade
  "Anemia",       140,    NA,     NA,     "G/L",  NA,       10,
  "Anemia",       139,    NA,     NA,     "G/L",  NA,       11,
  "Anemia",       100,    NA,     NA,     "G/L",  NA,       12,
  # ANRLO missing but AVAL satisfies grade 2
  "Anemia",       99,     NA,     NA,     "G/L",  "2",      13,
  # ANRLO missing but AVAL satisfies grade 3
  "Anemia",       79,     NA,     NA,     "G/L",  "3",      14,
  # Unit missing cannot grade
  "Anemia",       140,    140,    NA,     NA,     NA,       15,
  # AVAL missing cannot grade
  "Anemia",       NA,     140,    NA,     "G/L",  NA,       16,
)

input_anemia_si <- expected_anemia_si %>%
  select(-ATOXGRL)

## Test 5a: CTCAEv4 Anemia ----
test_that("derive_var_atoxgr_dir Test 5a: CTCAEv4 Anemia (SI unit)", {
  actual_anemia_ctcv4_si <- derive_var_atoxgr_dir(
    input_anemia_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_si,
    compare = actual_anemia_ctcv4_si,
    keys = c("TESTNUM")
  )
})

# CV_UNIT is g/dL so divide numeric values by 10
expected_anemia_cv <- expected_anemia_si %>%
  mutate(
    AVAL = if_else(is.na(AVAL), NA_real_, AVAL / 10),
    ANRLO = if_else(is.na(ANRLO), NA_real_, ANRLO / 10),
    ANRHI = if_else(is.na(ANRHI), NA_real_, ANRHI / 10),
    AVALU = if_else(toupper(AVALU) == "G/L", "g/dL", AVALU)
  )

input_anemia_cv <- expected_anemia_cv %>%
  select(-ATOXGRL)

## Test 5b: CTCAEv4 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 5b: CTCAEv4 Anemia (CV unit)", {
  actual_anemia_ctcv4 <- derive_var_atoxgr_dir(
    input_anemia_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_cv,
    compare = actual_anemia_ctcv4,
    keys = c("TESTNUM")
  )
})

## Test 5c: CTCAEv4 Anemia (SI + CV unit) ----
test_that("derive_var_atoxgr_dir Test 5c: CTCAEv4 Anemia (SI + CV unit)", {
  expected_anemia_sicv <- expected_anemia_cv %>%
    mutate(TESTNUM = TESTNUM + 16) %>%
    bind_rows(expected_anemia_si)

  input_anemia_sicv <- expected_anemia_sicv %>%
    select(-ATOXGRL)

  atoxgr_criteria_ctcv4_siuscv <- atoxgr_criteria_ctcv4_uscv %>%
    bind_rows(atoxgr_criteria_ctcv4)

  actual_anemia_ctcv4 <- derive_var_atoxgr_dir(
    input_anemia_sicv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_siuscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_sicv,
    compare = actual_anemia_ctcv4,
    keys = c("TESTNUM")
  )
})

## Test 6a: CTCAEv5 Anemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 6a: CTCAEv5 Anemia (SI unit)", {
  actual_anemia_ctcv5_si <- derive_var_atoxgr_dir(
    input_anemia_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_si,
    compare = actual_anemia_ctcv5_si,
    keys = c("TESTNUM")
  )
})

## Test 6b: CTCAEv5 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 6b: CTCAEv5 Anemia (CV unit)", {
  actual_anemia_ctcv5_cv <- derive_var_atoxgr_dir(
    input_anemia_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_cv,
    compare = actual_anemia_ctcv5_cv,
    keys = c("TESTNUM")
  )
})

## Test 6c: CTCAEv6 Anemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 6c: CTCAEv6 Anemia (SI unit)", {
  actual_anemia_ctcv6_si <- derive_var_atoxgr_dir(
    input_anemia_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_si,
    compare = actual_anemia_ctcv6_si,
    keys = c("TESTNUM")
  )
})

## Test 6d: CTCAEv6 Anemia (CV unit) ----
test_that("derive_var_atoxgr_dir Test 6d: CTCAEv6 Anemia (CV unit)", {
  actual_anemia_ctcv6_cv <- derive_var_atoxgr_dir(
    input_anemia_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia_cv,
    compare = actual_anemia_ctcv6_cv,
    keys = c("TESTNUM")
  )
})

### Leukocytosis
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)

### Grade 3: >100,000/mm3

expected_leukocytosis_si <- tibble::tribble(
  ~ATOXDSCH,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
  "Not a term",   99,     0,      NA,     "10^9/L",  NA,
  NA,             99,     0,      NA,     "10^9/L",  NA,
  "Leukocytosis", 101,    0,      40,     "10^9/L",  "3",
  "leukocytosis", 100,    0,      40,     "10^9/L",  "0",
  "Leukocytosis", 99,     0,      NA,     "10^9/L",  "0",
  # wrong UNIT - GRADE should be missing
  "Leukocytosis", 99,     0,      40,     "10^9/M",  NA,
  # Unit missing cannot grade
  "Leukocytosis", 99,     0,      40,     NA,        NA,
  # AVAL missing cannot grade
  "Leukocytosis", NA,     0,      40,     "10^9/L",  NA,
)
input_leukocytosis_si <- expected_leukocytosis_si %>%
  select(-ATOXGRH)


## Test 7a: CTCAEv4 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 7a: CTCAEv4 Leukocytosis (SI unit)", {
  actual_leukocytosis_si <- derive_var_atoxgr_dir(
    input_leukocytosis_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_si,
    compare = actual_leukocytosis_si,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

expected_leukocytosis_cv <- expected_leukocytosis_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

input_leukocytosis_cv <- expected_leukocytosis_cv %>%
  select(-ATOXGRH)

## Test 7b: CTCAEv4 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 7b: CTCAEv4 Leukocytosis (USCV unit)", {
  actual_leukocytosis_cv <- derive_var_atoxgr_dir(
    input_leukocytosis_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv,
    compare = actual_leukocytosis_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

expected_leukocytosis_cv2 <- expected_leukocytosis_si %>%
  mutate(
    AVAL = 1000 * AVAL,
    ANRHI = 1000 * ANRHI,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
input_leukocytosis_cv2 <- expected_leukocytosis_cv2 %>%
  select(-ATOXGRH)

## Test 7c: CTCAEv4 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 7c: CTCAEv4 Leukocytosis (legacy USCV unit)", {
  actual_leukocytosis_cv2 <- derive_var_atoxgr_dir(
    input_leukocytosis_cv2,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv2,
    compare = actual_leukocytosis_cv2,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8a: CTCAEv5 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 8a: CTCAEv5 Leukocytosis (SI unit)", {
  actual_leukocytosis_si <- derive_var_atoxgr_dir(
    input_leukocytosis_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_si,
    compare = actual_leukocytosis_si,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8b: CTCAEv5 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8b: CTCAEv5 Leukocytosis (USCV unit)", {
  actual_leukocytosis_cv <- derive_var_atoxgr_dir(
    input_leukocytosis_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv,
    compare = actual_leukocytosis_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8c: CTCAEv5 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8c: CTCAEv5 Leukocytosis (legacy USCV unit)", {
  actual_leukocytosis_cv2 <- derive_var_atoxgr_dir(
    input_leukocytosis_cv2,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv2,
    compare = actual_leukocytosis_cv2,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8d: CTCAEv6 Leukocytosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 8d: CTCAEv6 Leukocytosis (SI unit)", {
  actual_leukocytosisv6 <- derive_var_atoxgr_dir(
    input_leukocytosis_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_si,
    compare = actual_leukocytosisv6,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8e: CTCAEv6 Leukocytosis (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8e: CTCAEv6 Leukocytosis (USCV unit)", {
  actual_leukocytosis_cv <- derive_var_atoxgr_dir(
    input_leukocytosis_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv,
    compare = actual_leukocytosis_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 8f: CTCAEv6 Leukocytosis (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 8f: CTCAEv6 Leukocytosis (legacy USCV unit)", {
  actual_leukocytosis_cv2 <- derive_var_atoxgr_dir(
    input_leukocytosis_cv2,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis_cv2,
    compare = actual_leukocytosis_cv2,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})


## Investigations

### Activated partial thromboplastin time prolonged
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: >2.5 x ULN
### Grade 2: >1.5 - 2.5 x ULN
### Grade 1: >ULN - 1.5 x ULN

expected_aptt <- tibble::tribble(
  ~ATOXDSCH,                                         ~AVAL,  ~ANRHI,  ~AVALU,         ~ATOXGRH,
  "Not a term",                                      80,     100,     NA_character_,  NA,
  NA_character_,                                     60,     100,     NA_character_,  NA,
  "Activated partial thromboplastin time prolonged", 251,    100,     NA_character_,  "3",
  "Activated Partial thromboplastin time prolonged", 250,    100,     NA_character_,  "2",
  "Activated partial Thromboplastin time prolonged", 151,    100,     NA_character_,  "2",
  "Activated partial thromboplastin time prolonged", 150,    100,     NA_character_,  "1",
  "Activated partial thromboplastin Time prolonged", 101,    100,     NA_character_,  "1",
  "Activated partial thromboplastin time prolonged", 100,    100,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Activated partial thromboplastin time prolonged", 100,    NA,      NA_character_,  NA,
  # AVAL missing cannot grade
  "Activated partial thromboplastin time prolonged", NA,     100,     NA_character_,  NA,
)
input_aptt <- expected_aptt %>%
  select(-ATOXGRH)

## Test 9a: CTCAEv4 Activated partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 9: CTCAEv4 Activated partial thromboplastin time prolonged", {
  actual_aptt <- derive_var_atoxgr_dir(
    input_aptt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_aptt,
    compare = actual_aptt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 9b: CTCAEv4 Activated partial thromboplastin time prolonged (SI + CV METDATA) ----
test_that(paste(
  "derive_var_atoxgr_dir Test 9b: CTCAEv4 Activated partial thromboplastin time prolonged",
  "(SI + CV METDATA)"
), {
  atoxgr_criteria_ctcv4_sicv <- atoxgr_criteria_ctcv4 %>%
    bind_rows(atoxgr_criteria_ctcv4_uscv)

  actual_aptt <- derive_var_atoxgr_dir(
    input_aptt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_sicv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_aptt,
    compare = actual_aptt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 10a: CTCAEv5 Act. partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 10a: CTCAEv5 Act. partial thromboplastin time prolonged", {
  actual_aptt <- derive_var_atoxgr_dir(
    input_aptt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_aptt,
    compare = actual_aptt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 10b: CTCAEv6 Act. partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr_dir Test 10b: CTCAEv6 Act. partial thromboplastin time prolonged", {
  actual_aptt <- derive_var_atoxgr_dir(
    input_aptt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_aptt,
    compare = actual_aptt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 10c: CTCAEv6 Act. part. thromboplastin time prolonged (CV) ----
test_that("derive_var_atoxgr_dir Test 10c: CTCAEv6 Act. part. thromboplastin time prolonged (CV)", {
  actual_aptt <- derive_var_atoxgr_dir(
    input_aptt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_aptt,
    compare = actual_aptt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### Alanine aminotransferase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN

expected_alt_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                            ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",                         80,     40,     NA_character_,  NA,
  NA_character_,                        60,     40,     NA_character_,  NA,
  "Alanine aminotransferase Increased", 801,    40,     NA_character_,  "4",
  "Alanine aminotransferase Increased", 800,    40,     NA_character_,  "3",
  "Alanine aminotransferase Increased", 201,    40,     NA_character_,  "3",
  "Alanine aminotransferase Increased", 200,    40,     NA_character_,  "2",
  "Alanine aminotransferase Increased", 121,    40,     NA_character_,  "2",
  "Alanine aminotransferase Increased", 120,    40,     NA_character_,  "1",
  "Alanine aminotransferase Increased", 41,     40,     NA_character_,  "1",
  "Alanine aminotransferase Increased", 40,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Alanine aminotransferase Increased", 100,    NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Alanine aminotransferase Increased", NA,     40,     NA_character_,  NA,
)

## Test 11: CTCAEv4 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 11: CTCAEv4 Alanine aminotransferase increased", {
  input_alt <- expected_alt_ctcv4 %>%
    select(-ATOXGRH)

  actual_alt <- derive_var_atoxgr_dir(
    input_alt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alt_ctcv4,
    compare = actual_alt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})


### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >3.0 - 5.0 x ULN if BL was normal OR >3.0 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 3.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal

## Test 12a: CTCAEv5 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 12a: CTCAEv5 Alanine aminotransferase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_alt_ctcv5_norm <- expected_alt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_alt_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                            ~AVAL, ~BASE,        ~AVALU, ~ATOXGRH,
    "Not a term",                            80,    40, NA_character_,       NA,
    NA_character_,                           60,    40, NA_character_,       NA,
    "Alanine aminotransferase Increased",   801,    40, NA_character_,      "4",
    "Alanine aminotransferase Increased",   800,    40, NA_character_,      "3",
    "Alanine aminotransferase Increased",   201,    40, NA_character_,      "3",
    "Alanine aminotransferase Increased",   200,    40, NA_character_,      "2",
    "Alanine aminotransferase Increased",   121,    40, NA_character_,      "2",
    "Alanine aminotransferase Increased",   120,    40, NA_character_,      "1",
    "Alanine aminotransferase Increased",    60,    40, NA_character_,      "1",
    "Alanine aminotransferase Increased",    59,    40, NA_character_,      "0",
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased",   100,    NA, NA_character_,       NA,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased",    NA,    40, NA_character_,       NA,
  ) %>%
    # set BASE to be HIGH (not normal)
    mutate(
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  expected_alt_ctcv5 <- expected_alt_ctcv5_norm %>%
    bind_rows(expected_alt_ctcv5_abn)


  input_alt <- expected_alt_ctcv5 %>%
    select(-ATOXGRH)

  actual_alt <- derive_var_atoxgr_dir(
    input_alt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alt_ctcv5,
    compare = actual_alt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
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
  # V6 and V5 criteria identical when BASELINE normal
  expected_alt_ctcv6_norm <- expected_alt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_alt_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,                            ~AVAL, ~BASE,        ~AVALU, ~ATOXGRH,
    "Not a term",                            80,    40, NA_character_,       NA,
    NA_character_,                           60,    40, NA_character_,       NA,
    "Alanine aminotransferase Increased",   161,    40, NA_character_,      "4",
    "Alanine aminotransferase Increased",   160,    40, NA_character_,      "3",
    "Alanine aminotransferase Increased",    81,    40, NA_character_,      "3",
    "Alanine aminotransferase Increased",    80,    40, NA_character_,      "2",
    "Alanine aminotransferase Increased",    61,    40, NA_character_,      "2",
    "Alanine aminotransferase Increased",    60,    40, NA_character_,      "1",
    "Alanine aminotransferase Increased",    40,    40, NA_character_,      "1",
    "Alanine aminotransferase Increased",    39,    40, NA_character_,      "0",
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased",   100,    NA, NA_character_,       NA,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased",    NA,    40, NA_character_,       NA,
  ) %>%
    # set BASE to be HIGH (not normal)
    mutate(
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  expected_alt_ctcv6 <- expected_alt_ctcv6_norm %>%
    bind_rows(expected_alt_ctcv6_abn)


  input_alt <- expected_alt_ctcv6 %>%
    select(-ATOXGRH)

  actual_alt <- derive_var_atoxgr_dir(
    input_alt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alt_ctcv6,
    compare = actual_alt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
})

### Alkaline phosphatase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN
expected_alkph_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                         ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",                      80,     40,     NA_character_,  NA,
  NA_character_,                     60,     40,     NA_character_,  NA,
  "Alkaline phosphatase increased",  801,    40,     NA_character_,  "4",
  "Alkaline phosphatase increased",  800,    40,     NA_character_,  "3",
  "Alkaline phosphatase increased",  201,    40,     NA_character_,  "3",
  "Alkaline phosphatase increased",  200,    40,     NA_character_,  "2",
  "Alkaline phosphatase increased",  101,    40,     NA_character_,  "2",
  "Alkaline phosphatase increased",  100,    40,     NA_character_,  "1",
  "Alkaline phosphatase increased",  41,     40,     NA_character_,  "1",
  "Alkaline phosphatase increased",  40,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Alkaline phosphatase increased",  100,    NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Alkaline phosphatase increased",  NA,     40,     NA_character_,  NA,
)


## Test 13: CTCAEv4 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 13: CTCAEv4 Alkaline phosphatase increased", {
  input_alkph <- expected_alkph_ctcv4 %>%
    select(-ATOXGRH)

  actual_alkph <- derive_var_atoxgr_dir(
    input_alkph,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkph_ctcv4,
    compare = actual_alkph,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### Alkaline phosphatase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >2.5 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR >2.0 - 2.5 x BL if BL was abnormal

## Test 14a: CTCAEv5 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 14a: CTCAEv5 Alkaline phosphatase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_alkph_ctcv5_norm <- expected_alkph_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria

  expected_alkph_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                         ~AVAL,  ~BASE,  ~AVALU,         ~ATOXGRH,
    "Not a term",                      80,     40,     NA_character_,  NA,
    NA_character_,                     60,     40,     NA_character_,  NA,
    "Alkaline phosphatase increased",  801,    40,     NA_character_,  "4",
    "Alkaline phosphatase increased",  800,    40,     NA_character_,  "3",
    "Alkaline phosphatase increased",  201,    40,     NA_character_,  "3",
    "Alkaline phosphatase increased",  200,    40,     NA_character_,  "2",
    "Alkaline phosphatase increased",  101,    40,     NA_character_,  "2",
    "Alkaline phosphatase increased",  100,    40,     NA_character_,  "1",
    "Alkaline phosphatase increased",  80,     40,     NA_character_,  "1",
    "Alkaline phosphatase increased",  79,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Alkaline phosphatase increased",  100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Alkaline phosphatase increased",  NA,     40,     NA_character_,  NA,
  ) %>%
    # set BASE to be abnormal (HIGH) and create FLAG
    mutate(
      ANRHI = BASE - 1,
      BNRIND = "HIGH"
    )

  # combine records with baseline normal and abnormal
  expected_alkph_ctcv5 <- expected_alkph_ctcv5_norm %>%
    bind_rows(expected_alkph_ctcv5_abn)

  input_alkph <- expected_alkph_ctcv5 %>%
    select(-ATOXGRH)

  actual_alkph <- derive_var_atoxgr_dir(
    input_alkph,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkph_ctcv5,
    compare = actual_alkph,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
})

### Alkaline phosphatase increased
### NCICTCAEv6 different to NCICTCAEv5
### Grade 1: >Baseline and ULN

expected_alkph_ctcv6 <- tibble::tribble(
  ~ATOXDSCH,                         ~AVAL,  ~BASE,  ~ANRHI,  ~ATOXGRH,
  "Not a term",                      80,     40,     40,      NA,
  NA_character_,                     60,     40,     40,      NA,
  "Alkaline phosphatase increased",  41,     40,     40,      "1",
  "Alkaline phosphatase increased",  40,     40,     40,      "0",
  "Alkaline phosphatase increased",  41,     41,     40,      "0",
  "Alkaline phosphatase increased",  41,     40,     41,      "0",
  "Alkaline phosphatase increased",  40,     40,     NA,      "0",
  "Alkaline phosphatase increased",  40,     NA,     40,      "0",
  # ANRHI missing - cannot grade
  "Alkaline phosphatase increased",  100,    40,     NA,      NA,
  # BASE missing - cannot grade
  "Alkaline phosphatase increased",  41,     NA,     40,      NA,
  # AVAL missing cannot grade
  "Alkaline phosphatase increased",  NA,     40,     40,      NA,
) %>%
  mutate(AVALU = NA_character_)

input_alkph_ctcv6 <- expected_alkph_ctcv6 %>%
  select(-ATOXGRH)

## Test 14b: CTCAEv6 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr_dir Test 14b: CTCAEv6 Alkaline phosphatase increased", {
  actual_alkph_ctcv6 <- derive_var_atoxgr_dir(
    input_alkph_ctcv6,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkph_ctcv6,
    compare = actual_alkph_ctcv6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BASE")
  )
})

## Test 14c: CTCAEv6 Alkaline phosphatase increased (CV unit) ----
test_that("derive_var_atoxgr_dir Test 14b: CTCAEv6 Alkaline phosphatase increased (CV unit)", {
  actual_alkph_ctcv6 <- derive_var_atoxgr_dir(
    input_alkph_ctcv6,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkph_ctcv6,
    compare = actual_alkph_ctcv6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BASE")
  )
})


### Aspartate aminotransferase increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN

expected_ast_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                              ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",                           80,     40,     NA_character_,  NA,
  NA_character_,                          60,     40,     NA_character_,  NA,
  "Aspartate aminotransferase Increased", 801,    40,     NA_character_,  "4",
  "Aspartate aminotransferase Increased", 800,    40,     NA_character_,  "3",
  "Aspartate aminotransferase Increased", 201,    40,     NA_character_,  "3",
  "Aspartate aminotransferase Increased", 200,    40,     NA_character_,  "2",
  "Aspartate aminotransferase Increased", 121,    40,     NA_character_,  "2",
  "Aspartate aminotransferase Increased", 120,    40,     NA_character_,  "1",
  "Aspartate aminotransferase Increased", 41,     40,     NA_character_,  "1",
  "Aspartate aminotransferase Increased", 40,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Aspartate aminotransferase Increased", 100,    NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Aspartate aminotransferase Increased", NA,     40,     NA_character_,  NA,
)

## Test 15: CTCAEv4 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 15: CTCAEv4 Aspartate aminotransferase increased", {
  input_ast <- expected_ast_ctcv4 %>%
    select(-ATOXGRH)

  actual_ast <- derive_var_atoxgr_dir(
    input_ast,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ast_ctcv4,
    compare = actual_ast,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >3.0 - 5.0 x ULN if BL was normal OR >3.0 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 3.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal

## Test 16a: CTCAEv5 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr_dir Test 16a: CTCAEv5 Aspartate aminotransferase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_ast_ctcv5_norm <- expected_ast_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_ast_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,                              ~AVAL,  ~BASE,  ~AVALU,         ~ATOXGRH,
    "Not a term",                           80,     40,     NA_character_,  NA,
    NA_character_,                          60,     40,     NA_character_,  NA,
    "Aspartate aminotransferase Increased", 801,    40,     NA_character_,  "4",
    "Aspartate aminotransferase Increased", 800,    40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 201,    40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 200,    40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 121,    40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 120,    40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 60,     40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 59,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Aspartate aminotransferase Increased", 100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Aspartate aminotransferase Increased", NA,     40,     NA_character_,  NA,
  ) %>%
    # set BASE to be abnormal (HIGH or H)  and create FLAG
    mutate(
      ANRHI = BASE - 1,
      BNRIND = if_else(ATOXGRH %in% c("3", "4"), "H", "HIGH")
    )

  # combine records with baseline normal and abnormal
  expected_ast_ctcv5 <- expected_ast_ctcv5_norm %>%
    bind_rows(expected_ast_ctcv5_abn)

  input_ast <- expected_ast_ctcv5 %>%
    select(-ATOXGRH)

  actual_ast <- derive_var_atoxgr_dir(
    input_ast,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = c("HIGH", "H"),
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ast_ctcv5,
    compare = actual_ast,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
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
  # V5 and V4 criteria identical when BASELINE normal
  expected_ast_ctcv6_norm <- expected_ast_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_ast_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,                              ~AVAL,  ~BASE,  ~AVALU,         ~ATOXGRH,
    "Not a term",                           80,     40,     NA_character_,  NA,
    NA_character_,                          60,     40,     NA_character_,  NA,
    "Aspartate aminotransferase Increased", 161,    40,     NA_character_,  "4",
    "Aspartate aminotransferase Increased", 160,    40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 81,     40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 80,     40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 61,     40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 60,     40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 40,     40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 39,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Aspartate aminotransferase Increased", 100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Aspartate aminotransferase Increased", NA,     40,     NA_character_,  NA,
  ) %>%
    # set BASE to be abnormal (HIGH or H)  and create FLAG
    mutate(
      ANRHI = BASE - 1,
      BNRIND = if_else(ATOXGRH %in% c("3", "4"), "H", "HIGH")
    )

  # combine records with baseline normal and abnormal
  expected_ast_ctcv6 <- expected_ast_ctcv6_norm %>%
    bind_rows(expected_ast_ctcv6_abn)

  input_ast <- expected_ast_ctcv6 %>%
    select(-ATOXGRH)

  actual_ast <- derive_var_atoxgr_dir(
    input_ast,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = c("HIGH", "H"),
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ast_ctcv6,
    compare = actual_ast,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
})

### Blood bilirubin increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >10.0 x ULN
### Grade 3: >3.0 - 10.0 x ULN
### Grade 2: >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

expected_bili_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                   ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",                80,     40,     NA_character_,  NA,
  NA_character_,               60,     40,     NA_character_,  NA,
  "Blood bilirubin increased", 401,    40,     NA_character_,  "4",
  "Blood bilirubin increased", 400,    40,     NA_character_,  "3",
  "Blood bilirubin increased", 121,    40,     NA_character_,  "3",
  "Blood bilirubin increased", 120,    40,     NA_character_,  "2",
  "Blood bilirubin increased", 61,     40,     NA_character_,  "2",
  "Blood bilirubin increased", 60,     40,     NA_character_,  "1",
  "Blood bilirubin increased", 41,     40,     NA_character_,  "1",
  "Blood bilirubin increased", 40,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Blood bilirubin increased", 100,    NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Blood bilirubin increased", NA,     40,     NA_character_,  NA,
)

## Test 17: CTCAEv4 Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 17: CTCAEv4 Blood bilirubin increased", {
  input_bili <- expected_bili_ctcv4 %>%
    select(-ATOXGRH)

  actual_bili <- derive_var_atoxgr_dir(
    input_bili,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bili_ctcv4,
    compare = actual_bili,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 18a: CTCAEv5  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 18a: CTCAEv5  Blood bilirubin increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_bili_ctcv5_norm <- expected_bili_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with abnormal BASE then add records with normal BASE
  expected_bili_ctcv5 <- expected_bili_ctcv4 %>%
    # set BASE to ANRHI then make ANRHI < BASE
    mutate(
      BASE = ANRHI,
      ANRHI = ANRHI - 1,
      BNRIND = "HIGH"
    ) %>%
    bind_rows(expected_bili_ctcv5_norm)

  input_bili <- expected_bili_ctcv5 %>%
    select(-ATOXGRH)

  actual_bili <- derive_var_atoxgr_dir(
    input_bili,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bili_ctcv5,
    compare = actual_bili,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
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
expected_bili_ctcv6_norm <- expected_bili_ctcv4 %>%
  # set BASE to be normal (not HIGH) and create FLAG
  mutate(
    BASE = ANRHI,
    BNRIND = "NORMAL"
  )

# create records with abnormal BASE then add records with normal BASE
expected_bili_ctcv6 <- tibble::tribble(
  ~ATOXDSCH,                   ~AVAL,  ~BASE,  ~AVALU,         ~ATOXGRH,
  "Not a term",                80,     40,     NA_character_,  NA,
  NA_character_,               60,     40,     NA_character_,  NA,
  "Blood bilirubin increased", 401,    40,     NA_character_,  "4",
  "Blood bilirubin increased", 400,    40,     NA_character_,  "3",
  "Blood bilirubin increased", 101,    40,     NA_character_,  "3",
  "Blood bilirubin increased", 100,    40,     NA_character_,  "2",
  "Blood bilirubin increased", 61,     40,     NA_character_,  "2",
  "Blood bilirubin increased", 60,     40,     NA_character_,  "1",
  "Blood bilirubin increased", 40,     40,     NA_character_,  "1",
  "Blood bilirubin increased", 39,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Blood bilirubin increased", 100,    NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Blood bilirubin increased", NA,     40,     NA_character_,  NA,
) %>%
  # set BASE to ANRHI then make ANRHI < BASE
  mutate(
    ANRHI = BASE - 1,
    BNRIND = "HIGH"
  ) %>%
  bind_rows(expected_bili_ctcv6_norm)

input_bili_ctcv6 <- expected_bili_ctcv6 %>%
  select(-ATOXGRH)

## Test 18b: CTCAEv6  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 18b: CTCAEv6  Blood bilirubin increased", {
  actual_bili_ctcv6 <- derive_var_atoxgr_dir(
    input_bili_ctcv6,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bili_ctcv6,
    compare = actual_bili_ctcv6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
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

expected_cd4_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~AVALU,    ~ATOXGRL,
  "Not a term",                 80,     120,    "10^9/L",  NA,
  NA_character_,                60,     50,     "10^9/L",  NA,
  "CD4 lymphocytes decreased",  0.04,   0.8,    "10^9/L",  "4",
  "CD4 lymphocytes decreased",  0.05,   0.8,    "10^9/l",  "3",
  "CD4 lymphocytes decreased",  0.19,   0.8,    "10^9/L",  "3",
  "CD4 lymphocytes decreased",  0.2,    0.8,    "10^9/L",  "2",
  "CD4 lymphocytes decreased",  0.49,   0.8,    "10^9/L",  "2",
  # wrong unit - grade missing
  "CD4 lymphocytes decreased",  0.49,   0.8,    "10^8/L",  NA,
  "CD4 lymphocytes decreased",  0.5,    0.8,    "10^9/L",  "1",
  "CD4 lymphocytes decreased",  0.79,   0.8,    "10^9/L",  "1",
  "CD4 lymphocytes decreased",  0.8,    0.8,    "10^9/L",  "0",
  # ANRLO missing - AVAL satisfies grade 2 - 4
  "CD4 lymphocytes decreased",  0.04,   NA,     "10^9/L",  "4",
  "CD4 lymphocytes decreased",  0.05,   NA,     "10^9/L",  "3",
  "CD4 lymphocytes decreased",  0.19,   NA,     "10^9/L",  "3",
  "CD4 lymphocytes decreased",  0.2,    NA,     "10^9/L",  "2",
  "CD4 lymphocytes decreased",  0.49,   NA,     "10^9/L",  "2",
  # ANRLO missing - AVAL does NOT satisfies grade 2 - 4
  "CD4 lymphocytes decreased",  0.5,    NA,     "10^9/L",  NA,
  "CD4 lymphocytes decreased",  0.79,   NA,     "10^9/L",  NA,
  "CD4 lymphocytes decreased",  0.8,    NA,     "10^9/L",  NA,
  # Unit missing - cannot grade
  "CD4 lymphocytes decreased",  0.8,    0.8,    NA,        NA,
  # AVAL missing cannot grade
  "CD4 lymphocytes decreased",  NA,     0.8,    "10^9/L",  NA,
)

input_cd4_si <- expected_cd4_si %>%
  select(-ATOXGRL)

expected_cd4_cv <- expected_cd4_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "1/uL", NA_character_)
  )

input_cd4_cv <- expected_cd4_cv %>%
  select(-ATOXGRL)

## Test 19a: CTCAEv4 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 19a: CTCAEv4 CD4 Lymphocytes decreased (SI unit)", {
  actual_cd4_si <- derive_var_atoxgr_dir(
    input_cd4_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_si,
    compare = actual_cd4_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

## Test 19b: CTCAEv4 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 19b: CTCAEv4 CD4 Lymphocytes decreased (USCV unit)", {
  actual_cd4_cv <- derive_var_atoxgr_dir(
    input_cd4_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_cv,
    compare = actual_cd4_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

## Test 20a: CTCAEv5 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 20a: CTCAEv5 CD4 Lymphocytes decreased (SI unit)", {
  actual_cd4_si <- derive_var_atoxgr_dir(
    input_cd4_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_si,
    compare = actual_cd4_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

## Test 20b: CTCAEv5 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 20b: CTCAEv5 CD4 Lymphocytes decreased (USCV unit)", {
  actual_cd4_cv <- derive_var_atoxgr_dir(
    input_cd4_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_cv,
    compare = actual_cd4_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

## Test 20c: CTCAEv6 CD4 Lymphocytes decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 20c: CTCAEv6 CD4 Lymphocytes decreased (SI unit)", {
  actual_cd4_si <- derive_var_atoxgr_dir(
    input_cd4_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_si,
    compare = actual_cd4_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

## Test 20d: CTCAEv6 CD4 Lymphocytes decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 20d: CTCAEv6 CD4 Lymphocytes decreased (USCV unit)", {
  actual_cd4_cv <- derive_var_atoxgr_dir(
    input_cd4_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4_cv,
    compare = actual_cd4_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
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

expected_choles_si <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
  "Not a term",       8,      0,      5,      "mmol/L",  NA,
  NA_character_,      10,     0,      5,      "mmol/L",  NA,
  "Cholesterol high", 12.93,  0,      5,      "mmol/L",  "4",
  "Cholesterol High", 12.92,  0,      5,      "mmol/L",  "3",
  "Cholesterol high", 10.35,  0,      5,      "Mmol/L",  "3",
  # wrong unit - grade missing
  "Cholesterol high", 10.35,  0,      5,      "umol/L",  NA,
  "Cholesterol high", 10.34,  0,      5,      "mmol/L",  "2",
  "Cholesterol high", 7.76,   0,      5,      "mmol/L",  "2",
  "Cholesterol high", 7.75,   0,      5,      "mmol/L",  "1",
  "Cholesterol high", 5.1,    0,      5,      "mmol/L",  "1",
  "Cholesterol high", 5,      0,      5,      "mmol/L",  "0",
  # ANRHI missing - AVAL satisfies grade 2 - 4
  "Cholesterol high", 12.93,  0,      NA,     "mmol/L",  "4",
  "Cholesterol High", 12.92,  0,      NA,     "mmol/L",  "3",
  "Cholesterol high", 10.35,  0,      NA,     "Mmol/L",  "3",
  "Cholesterol high", 10.34,  0,      NA,     "mmol/L",  "2",
  "Cholesterol high", 7.76,   0,      NA,     "mmol/L",  "2",
  # ANRHI missing - AVAL does NOT satisfies grade 2 - 4
  "Cholesterol high", 7.75,   0,      NA,     "mmol/L",  NA,
  "Cholesterol high", 5.1,    0,      NA,     "mmol/L",  NA,
  "Cholesterol high", 5,      0,      NA,     "mmol/L",  NA,
  # Unit missing - cannot grade
  "Cholesterol high", 5,      0,      5,      NA,        NA,
  # AVAL missing cannot grade
  "Cholesterol high", NA,     0,      5,      "mmol/L",  NA,
)

expected_choles_cv <- tibble::tribble(
  ~ATOXDSCH,           ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH,
  "Not a term",        8,      0,       200,     "mg/dL",        NA,
  NA_character_,       10,     0,       200,     "mg/dL",        NA,
  "Cholesterol high",  501,    0,       200,     "mg/dL",       "4",
  "Cholesterol High",  500,    0,       200,     "mg/dL",       "3",
  "Cholesterol high",  401,    0,       200,     "Mg/dL",       "3",
  # wrong unit - grade missing
  "Cholesterol high",  400,    0,       200,     "ug/dL",        NA,
  "Cholesterol high",  400,    0,       200,     "mg/dL",       "2",
  "Cholesterol high",  301,    0,       200,     "mg/dL",       "2",
  "Cholesterol high",  300,    0,       200,     "mg/dL",       "1",
  "Cholesterol high",  201,    0,       200,     "mg/dL",       "1",
  "Cholesterol high",  200,    0,       200,     "mg/dL",       "0",
  # ANRHI missing - AVAL satisfies grade 2 - 4
  "Cholesterol high",  501,    0,       NA,      "mg/dL",       "4",
  "Cholesterol High",  500,    0,       NA,      "mg/dL",       "3",
  "Cholesterol high",  401,    0,       NA,      "MG/dL",       "3",
  "Cholesterol high",  400,    0,       NA,      "mg/dL",       "2",
  "Cholesterol high",  301,    0,       NA,      "mg/dL",       "2",
  # ANRHI missing - AVAL does NOT satisfies grade 2 - 4
  "Cholesterol high",  300,    0,       NA,      "mg/dL",        NA,
  "Cholesterol high",  201,    0,       NA,      "mg/dL",        NA,
  "Cholesterol high",  200,    0,       NA,      "mg/dL",        NA,
  # Unit missing - cannot grade
  "Cholesterol high",  200,    0,       200,          NA,        NA,
  # AVAL missing cannot grade
  "Cholesterol high",  NA,     0,       200,    "mmol/L",        NA,
)


input_choles_si <- expected_choles_si %>%
  select(-ATOXGRH)

input_choles_cv <- expected_choles_cv %>%
  select(-ATOXGRH)

## Test 21a: CTCAEv4 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 21a: CTCAEv4 Cholesterol high (SI unit)", {
  actual_choles_si <- derive_var_atoxgr_dir(
    input_choles_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_si,
    compare = actual_choles_si,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 21b: CTCAEv4 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 21b: CTCAEv4 Cholesterol high (USCV unit)", {
  actual_choles_cv <- derive_var_atoxgr_dir(
    input_choles_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_cv,
    compare = actual_choles_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 22a: CTCAEv5 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 22a: CTCAEv5 Cholesterol high (SI unit)", {
  actual_choles_si <- derive_var_atoxgr_dir(
    input_choles_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_si,
    compare = actual_choles_si,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 22b: CTCAEv5 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 22b: CTCAEv5 Cholesterol high (USCV unit)", {
  actual_choles_cv <- derive_var_atoxgr_dir(
    input_choles_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_cv,
    compare = actual_choles_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 22c: CTCAEv6 Cholesterol high (SI unit) ----
test_that("derive_var_atoxgr_dir Test 22c: CTCAEv6 Cholesterol high (SI unit)", {
  actual_choles_si <- derive_var_atoxgr_dir(
    input_choles_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_si,
    compare = actual_choles_si,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 22d: CTCAEv6 Cholesterol high (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 22d: CTCAEv6 Cholesterol high (USCV unit)", {
  actual_choles_cv <- derive_var_atoxgr_dir(
    input_choles_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles_cv,
    compare = actual_choles_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### CPK increased
### NCICTCAEv5 same criteria as NCICTCAEv4
### Grade 4: >10.0 x ULN
### Grade 3: >5.0 - 10.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN

expected_cpk <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,  ~AVALU,         ~ATOXGRH,
  "Not a term",     80,     0,       40,      NA_character_,  NA,
  NA_character_,    60,     0,       40,      NA_character_,  NA,
  "CPK increased",  401,    0,       40,      NA_character_,  "4",
  "CPK increased",  400,    0,       40,      NA_character_,  "3",
  "CPK increased",  201,    0,       40,      NA_character_,  "3",
  "CPK increased",  200,    0,       40,      NA_character_,  "2",
  "CPK increased",  101,    0,       40,      NA_character_,  "2",
  "CPK increased",  100,    0,       40,      NA_character_,  "1",
  "CPK increased",  41,     0,       40,      NA_character_,  "1",
  "CPK increased",  40,     0,       40,      NA_character_,  "0",
  # ANRHI missing - cannot grade
  "CPK increased",  100,    0,       NA,      NA_character_,  NA,
  # AVAL missing cannot grade
  "CPK increased",  NA,     0,       40,      NA_character_,  NA,
)
input_cpk <- expected_cpk %>%
  select(-ATOXGRH)

## Test 23: CTCAEv4 CPK increased ----
test_that("derive_var_atoxgr_dir Test 23: CTCAEv4 CPK increased", {
  actual_cpk <- derive_var_atoxgr_dir(
    input_cpk,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cpk,
    compare = actual_cpk,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 24: CTCAEv5 CPK increased ----
test_that("derive_var_atoxgr_dir Test 24: CTCAEv5 CPK increased", {
  actual_cpk <- derive_var_atoxgr_dir(
    input_cpk,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cpk,
    compare = actual_cpk,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Creatinine increased (NCICTCv4)
### NCICTCAEv5 same criteria as NCICTCAEv4 except for Grade 1
### Grade 4: >6.0 x ULN
### Grade 3: >3.0 baseline; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x baseline; >1.5 - 3.0 x ULN
### Grade 1: >1 - 1.5 x baseline; >ULN - 1.5 x ULN

# create flag to remove obs not relevant for NCI-CTCAEv5
expected_creatn <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~ATOXGRH, ~V4, ~V5,
  "Not a term",            80,     80,    40,     NA_character_,  NA,       "Y", "Y",
  NA_character_,           60,     60,    40,     NA_character_,  NA,       "Y", "Y",
  # GRADE derived from AVAL against ANRHI
  "Creatinine increased",  241,    241,   40,     NA_character_,  "4",      "Y", "Y",
  "Creatinine increased",  240,    230,   40,     NA_character_,  "3",      "Y", "Y",
  "Creatinine increased",  121,    120,   40,     NA_character_,  "3",      "Y", "Y",
  "Creatinine increased",  120,    119,   40,     NA_character_,  "2",      "Y", "Y",
  "Creatinine increased",  61,     60,    40,     NA_character_,  "2",      "Y", "Y",
  "Creatinine increased",  60,     60,    40,     NA_character_,  "1",      "Y", "Y",
  "Creatinine increased",  41,     41,    40,     NA_character_,  "1",      "Y", "Y",
  "Creatinine increased",  40,     40,    40,     NA_character_,  "0",      "Y", "Y",
  # GRADE derived from AVAL against BASE
  "Creatinine increased",  42,     6,     40,     NA_character_,  "3",      "Y", "Y",
  "Creatinine increased",  42,     13.9,  40,     NA_character_,  "3",      "Y", "Y",
  "Creatinine increased",  42,     14,    40,     NA_character_,  "2",      "Y", "Y",
  "Creatinine increased",  42.1,   28,    40,     NA_character_,  "2",      "Y", "Y",
  "Creatinine increased",  42,     28,    42,     NA_character_,  "1",      "Y", "N",
  "Creatinine increased",  42,     41,    42,     NA_character_,  "1",      "Y", "N",
  "Creatinine increased",  42,     42,    42,     NA_character_,  "0",      "Y", "N",
  # BASE missing - AVAL <= ANRLO cannot grade as NORMAL
  "Creatinine increased",  42,     NA,    42,     NA_character_,  NA,       "Y", "N",
  # ANRHI missing - AVAL <= BASE cannot grade as NORMAL
  "Creatinine increased",  42,     42,    NA,     NA_character_,  NA,       "Y", "Y",
  # AVAL missing cannot grade
  "Creatinine increased",  NA,     0,     40,     NA_character_,  NA,       "Y", "Y",
)

## Test 25: CTCAEv4 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 25: CTCAEv4 Creatinine increased", {
  input_creatn <- expected_creatn %>%
    select(-ATOXGRH)

  actual_creatn <- derive_var_atoxgr_dir(
    input_creatn,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_creatn,
    compare = actual_creatn,
    keys = c("ATOXDSCH", "AVAL", "BASE", "ANRHI", "AVALU")
  )
})

### Creatinine increased (NCICTCv5)
### NCICTCAEv5 same criteria as NCICTCAEv4 except for Grade 1
### Grade 4: >6.0 x ULN
### Grade 3: >3.0 baseline; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x baseline; >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

## Test 26a: CTCAEv5 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 26a: CTCAEv5 Creatinine increased", {
  expected_creatn <- expected_creatn %>%
    filter(V5 == "Y")

  input_creatn <- expected_creatn %>%
    select(-ATOXGRH)

  actual_creatn <- derive_var_atoxgr_dir(
    input_creatn,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_creatn,
    compare = actual_creatn,
    keys = c("ATOXDSCH", "AVAL", "BASE", "ANRHI", "AVALU")
  )
})

### Creatinine increased (NCICTCv6)

### Grade 4: >6.0 x ULN
### Grade 3: >3.0 x if baseline is below LLN; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x if baseline is below LLN; >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

## Test 26b: CTCAEv6 Creatinine increased ----
test_that("derive_var_atoxgr_dir Test 26b: CTCAEv6 Creatinine increased", {
  expected_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~BNRIND,        ~ATOXGRH,
    "Not a term",            80,     80,    40,     NA_character_,  "NORMAL",       NA_character_,
    NA_character_,           60,     60,    40,     NA_character_,  "NORMAL",       NA_character_,
    # GRADE derived from AVAL against ANRHI
    "Creatinine increased",  241,    40,    40,     NA_character_,  "NORMAL",       "4",
    "Creatinine increased",  240,    40,    40,     NA_character_,  "NORMAL",       "3",
    "Creatinine increased",  121,    40,    40,     NA_character_,  "NORMAL",       "3",
    "Creatinine increased",  120,    40,    40,     NA_character_,  "NORMAL",       "2",
    "Creatinine increased",  61,     40,    40,     NA_character_,  "NORMAL",       "2",
    "Creatinine increased",  60,     40,    40,     NA_character_,  "NORMAL",       "1",
    "Creatinine increased",  41,     40,    40,     NA_character_,  "NORMAL",       "1",
    "Creatinine increased",  40,     40,    40,     NA_character_,  "NORMAL",       "0",
    # GRADE derived from AVAL against BASE when baseline is below LLN
    "Creatinine increased",  121,    40,    41,     NA_character_,  "LOW",          "3",
    "Creatinine increased",  120,    40,    41,     NA_character_,  "LOW",          "2",
    "Creatinine increased",  61,     40,    41,     NA_character_,  "LOW",          "2",
    "Creatinine increased",  60,     40,    60,     NA_character_,  "LOW",          "0",
    # BASE missing - AVAL <= ANRLO can grade as NORMAL
    "Creatinine increased",  42,     NA,    42,     NA_character_,  NA_character_,  "0",
    # ANRHI missing - baseline is below LLN and AVAL <= 1.5 x BASE but cannot grade as NORMAL
    "Creatinine increased",  60,     40,    NA,     NA_character_,  "LOW",          NA_character_,
    # AVAL missing cannot grade
    "Creatinine increased",  NA,     0,     40,     NA_character_,  "NORMAL",       NA_character_,
  )

  input_creatn <- expected_creatn %>%
    select(-ATOXGRH)

  actual_creatn <- derive_var_atoxgr_dir(
    input_creatn,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH",
    low_indicator = "LOW",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_creatn,
    compare = actual_creatn,
    keys = c("ATOXDSCH", "AVAL", "BASE", "ANRHI", "BNRIND")
  )
})

### Fibrinogen decreased (NCICTCv4)
### Grade 4: <0.25 x LLN or 75% decrease from baseline or absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN or 50 - <75% decrease from baseline
### Grade 2: <0.75 - 0.5 x LLN or 25 - <50% decrease from baseline
### Grade 1: <1.0 - 0.75 x LLN or <25% decrease from baseline

## Test 27: CTCAEv4 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 27: CTCAEv4 Fibrinogen decreased", {
  expected_fib <- tibble::tribble(
    ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~PCHG,  ~AVALU,  ~ATOXGRL,
    "Not a term",            9,      10,     40,     "g/L",   NA,
    NA_character_,           10,     10,     40,     "g/L",   NA,
    # Satisfies < 0.5 for grade 4 - other criteria missing
    "Fibrinogen decreased",  0.49,   NA,     NA,     "g/L",   "4",
    # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
    "Fibrinogen decreased",  0.49,   1,      -51,    "g/L",   "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG missing
    "Fibrinogen decreased",  0.5,    2.1,    NA,     "g/L",   "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
    "Fibrinogen decreased",  0.5,    2.1,    -51,    "g/L",   "4",
    # Satisfies <=75% decrease for grade 4 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -75,    "g/L",   "4",
    # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
    "Fibrinogen decreased",  1,      0.49,   -75,    "g/L",   "4",
    # Satisfies < 0.5*LLN for grade 3 - PCHG missing
    "Fibrinogen decreased",  1,      2.1,    NA,     "g/L",   "3",
    # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
    "Fibrinogen decreased",  1,      2.1,    -49,    "g/L",   "3",
    # Satisfies <=50% decrease for grade 3 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -50,    "g/L",   "3",
    # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
    "Fibrinogen decreased",  1,      2,      -50,    "g/L",   "3",
    # Satisfies < 0.75*LLN for grade 2 - PCHG missing
    "Fibrinogen decreased",  1.5,    2.1,    NA,     "g/L",   "2",
    # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
    "Fibrinogen decreased",  1.5,    2.1,    -10,    "g/L",   "2",
    # Satisfies <=25% for grade 2 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -25,    "g/L",   "2",
    # Satisfies <=25% for grade 2 - LLN satisfies grade 1
    "Fibrinogen decreased",  1.5,    1.6,    -25,    "g/L",   "2",
    # Satisfies < LLN for grade 1 - PCHG missing
    "Fibrinogen decreased",  2,      2.1,    NA,     "g/L",   "1",
    # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
    "Fibrinogen decreased",  2,      2.1,    10,     "g/L",   "1",
    # Satisfies % decrease for grade 1 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -1,     "g/L",   "1",
    # Satisfies % decrease for grade 1 - AVAL = LLN
    "Fibrinogen decreased",  1.5,    1.5,    -1,     "g/L",   "1",
    # Satisfies grade 0 - AVAL >= LLN AND no % descrease
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/L",   "0",
    # AVAL >= LLN BUT PCT missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    1.5,    NA,     "g/L",   NA,
    # PCT >= 0 BUT LLN missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    NA,     10,     "g/L",   NA,
    # AVAL missing cannot grade
    "Fibrinogen decreased",  NA,     1.5,    10,     "g/L",   NA,
    # wrong unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/dL",  NA,
    # missing unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      NA,      NA,
  )
  input_fib <- expected_fib %>%
    select(-ATOXGRL)

  actual_fib <- derive_var_atoxgr_dir(
    input_fib,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_fib,
    compare = actual_fib,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "PCHG", "AVALU")
  )
})

### Fibrinogen decreased (NCICTCv5 + NCICTCv6)
### Grade 4: <0.25 x LLN OR if abnormal, 75% dec. from BL OR absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN OR if abnormal, 50 - <75% dec. from BL
### Grade 2: <0.75 - 0.5 x LLN OR if abnormal, 25 - <50% dec. from BL
### Grade 1: <1.0 - 0.75 x LLN OR if abnormal, <25% dec. from BL

expected_fib <- tibble::tribble(
  ~ATOXDSCL,               ~AVAL,  ~ANRLO,  ~PCHG,  ~AVALU,  ~ATOXGRL,
  "Not a term",            9,      10,      40,      "g/L",        NA,
  NA_character_,           10,     10,      40,      "g/L",        NA,
  # Satisfies < 0.5 for grade 4 - other criteria missing
  "Fibrinogen decreased",  0.49,   NA,      NA,      "g/L",       "4",
  # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
  "Fibrinogen decreased",  0.49,   1,       -51,     "g/L",       "4",
  # Satisfies < 0.25*LLN for grade 4 - PCHG missing
  "Fibrinogen decreased",  0.5,    2.1,     NA,      "g/L",       "4",
  # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
  "Fibrinogen decreased",  0.5,    2.1,     -51,     "g/L",       "4",
  # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
  "Fibrinogen decreased",  0.5,    1.1,     -75,     "g/L",       "4",
  # Satisfies < 0.5*LLN for grade 3 - PCHG missing
  "Fibrinogen decreased",  1,      2.1,     NA,      "g/L",       "3",
  # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
  "Fibrinogen decreased",  1,      2.1,     -49,     "g/L",       "3",
  # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
  "Fibrinogen decreased",  1,      2,       -50,     "g/L",       "3",
  # Satisfies < 0.75*LLN for grade 2 - PCHG missing
  "Fibrinogen decreased",  1.5,    2.1,     NA,      "g/L",       "2",
  # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
  "Fibrinogen decreased",  1.5,    2.1,     -10,     "g/L",       "2",
  # Satisfies <=25% for grade 2 - LLN satisfies grade 1
  "Fibrinogen decreased",  1.5,    1.6,     -25,     "g/L",       "2",
  # Satisfies < LLN for grade 1 - PCHG missing
  "Fibrinogen decreased",  2,      2.1,     NA,      "g/L",       "1",
  # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
  "Fibrinogen decreased",  2,      2.1,     10,      "g/L",       "1",
  # Satisfies grade 0 - AVAL >= LLN AND no % descrease
  "Fibrinogen decreased",  1.5,    1.5,     0,       "g/L",       "0",
  # Satisfies % decrease for grade 1 - AVAL = LLN so not abnormal
  "Fibrinogen decreased",  1.5,    1.5,     -1,      "g/L",       "0",
  # AVAL >= LLN - PCT missing but its normal so ignore PCT
  "Fibrinogen decreased",  1.5,    1.5,     NA,      "g/L",       "0",
  # Satisfies <=75% decrease for grade 4 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1,      NA,      -75,     "g/L",        NA,
  # Satisfies <=50% decrease for grade 3 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1,      NA,      -50,     "g/L",        NA,
  # Satisfies <=25% decrease for grade 2 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1.5,    NA,      -25,     "g/L",        NA,
  # Satisfies % decrease for grade 1 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1.5,    NA,      -1,      "g/L",        NA,
  # PCT >= 0 BUT LLN missing cannot grade as NORMAL
  "Fibrinogen decreased",  1.5,    NA,      10,      "g/L",        NA,
  # AVAL missing cannot grade
  "Fibrinogen decreased",  NA,     1.5,     10,      "g/L",        NA,
  # wrong unit cannot grade as it may satisfy grade 4
  "Fibrinogen decreased",  1.5,    1.5,     0,      "g/dL",        NA,
  # missing unit cannot grade as it may satisfy grade 4
  "Fibrinogen decreased",  1.5,    1.5,     0,          NA,        NA,
)
input_fib <- expected_fib %>%
  select(-ATOXGRL)

## Test 28a: CTCAEv5 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 28a: CTCAEv5 Fibrinogen decreased", {
  actual_fib <- derive_var_atoxgr_dir(
    input_fib,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_fib,
    compare = actual_fib,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "PCHG", "AVALU")
  )
})

## Test 28b: CTCAEv6 Fibrinogen decreased ----
test_that("derive_var_atoxgr_dir Test 28b: CTCAEv6 Fibrinogen decreased", {
  actual_fib <- derive_var_atoxgr_dir(
    input_fib,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_fib,
    compare = actual_fib,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "PCHG", "AVALU")
  )
})

### GGT increased (NCICTCv4)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN

expected_ggt_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,        ~ATOXGRH,
  "Not a term",    80,    0,      40,     NA_character_, NA,
  NA_character_,   60,    0,      40,     NA_character_, NA,
  "GGT increased", 801,   0,      40,     NA_character_, "4",
  "GGT increased", 800,   0,      40,     NA_character_, "3",
  "GGT increased", 201,   0,      40,     NA_character_, "3",
  "GGT increased", 200,   0,      40,     NA_character_, "2",
  "GGT increased", 101,   0,      40,     NA_character_, "2",
  "GGT increased", 100,   0,      40,     NA_character_, "1",
  "GGT increased", 41,    0,      40,     NA_character_, "1",
  "GGT increased", 40,    0,      40,     NA_character_, "0",
  # ANRHI missing - cannot grade
  "GGT increased", 100,   0,      NA,     NA_character_, NA,
  # AVAL missing cannot grade
  "GGT increased", NA,    0,      NA,     NA_character_, NA,
)

## Test 29: CTCAEv4 GGT increased ----
test_that("derive_var_atoxgr_dir Test 29: CTCAEv4 GGT increased", {
  input_ggt <- expected_ggt_ctcv4 %>%
    select(-ATOXGRH)

  actual_ggt <- derive_var_atoxgr_dir(
    input_ggt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ggt_ctcv4,
    compare = actual_ggt,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### GGT increased (NCICTCv5)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >20.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >5.0 - 20.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >2.5 - 5.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR >2.0 - 2.5 x BL if BL was abnormal

## Test 30a: CTCAEv5 GGT increased ----
test_that("derive_var_atoxgr_dir Test 30a: CTCAEv5 GGT increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_ggt_ctcv5_norm <- expected_ggt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_ggt_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRLO, ~BASE,  ~AVALU,        ~ATOXGRH,
    "Not a term",    80,    0,      40,     NA_character_, NA,
    NA_character_,   60,    0,      40,     NA_character_, NA,
    "GGT increased", 801,   0,      40,     NA_character_, "4",
    "GGT increased", 800,   0,      40,     NA_character_, "3",
    "GGT increased", 201,   0,      40,     NA_character_, "3",
    "GGT increased", 200,   0,      40,     NA_character_, "2",
    "GGT increased", 101,   0,      40,     NA_character_, "2",
    "GGT increased", 100,   0,      40,     NA_character_, "1",
    "GGT increased", 81,    0,      40,     NA_character_, "1",
    "GGT increased", 80,    0,      40,     NA_character_, "0",
    # ANRHI missing - cannot grade
    "GGT increased", 100,   0,      NA,     NA_character_, NA,
    # AVAL missing cannot grade
    "GGT increased", NA,    0,      NA,     NA_character_, NA,
  ) %>%
    # set BASE to be abnormal (HIGH HIGH) and create FLAG
    # set high_indicator to "HIGH HIGH" also in function call
    mutate(
      ANRHI = BASE - 1,
      BNRIND = "HIGH HIGH"
    )

  # combine records with baseline normal and abnormal
  expected_ggt_ctcv5 <- expected_ggt_ctcv5_norm %>%
    bind_rows(expected_ggt_ctcv5_abn)

  input_ggt <- expected_ggt_ctcv5 %>%
    select(-ATOXGRH)

  actual_ggt <- derive_var_atoxgr_dir(
    input_ggt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ggt_ctcv5,
    compare = actual_ggt,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "BNRIND", "AVALU")
  )
})

### GGT increased (NCICTCv6)
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >20.0 x ULN if BL was normal OR >10.0 x BL if BL was abnormal
### Grade 3: >5.0 - 20.0 x ULN if BL was normal OR >3.0 - 10.0 x BL if BL was abnormal
### Grade 2: >2.5 - 5.0 x ULN if BL was normal OR >1.5 - 3.0 x BL if BL was abnormal
### Grade 1: >ULN - 2.5 x ULN if BL was normal OR 1.0 - 1.5 x BL if BL was abnormal

## Test 30b: CTCAEv6 GGT increased ----
test_that("derive_var_atoxgr_dir Test 30b: CTCAEv6 GGT increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_ggt_ctcv6_norm <- expected_ggt_ctcv4 %>%
    # set BASE to be normal (not HIGH) and create FLAG
    mutate(
      BASE = ANRHI,
      BNRIND = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_ggt_ctcv6_abn <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRLO, ~BASE,  ~AVALU,        ~ATOXGRH,
    "Not a term",    80,    0,      40,     NA_character_, NA,
    NA_character_,   60,    0,      40,     NA_character_, NA,
    "GGT increased", 401,   0,      40,     NA_character_, "4",
    "GGT increased", 400,   0,      40,     NA_character_, "3",
    "GGT increased", 121,   0,      40,     NA_character_, "3",
    "GGT increased", 120,   0,      40,     NA_character_, "2",
    "GGT increased", 61,    0,      40,     NA_character_, "2",
    "GGT increased", 60,    0,      40,     NA_character_, "1",
    "GGT increased", 40,    0,      40,     NA_character_, "1",
    "GGT increased", 39,    0,      40,     NA_character_, "0",
    # ANRHI missing - cannot grade
    "GGT increased", 100,   0,      NA,     NA_character_, NA,
    # AVAL missing cannot grade
    "GGT increased", NA,    0,      NA,     NA_character_, NA,
  ) %>%
    # set BASE to be abnormal (HIGH HIGH) and create FLAG
    # set high_indicator to "HIGH HIGH" also in function call
    mutate(
      ANRHI = BASE - 1,
      BNRIND = "HIGH HIGH"
    )

  # combine records with baseline normal and abnormal
  expected_ggt_ctcv6 <- expected_ggt_ctcv6_norm %>%
    bind_rows(expected_ggt_ctcv6_abn)

  input_ggt <- expected_ggt_ctcv6 %>%
    select(-ATOXGRH)

  actual_ggt <- derive_var_atoxgr_dir(
    input_ggt,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    high_indicator = "HIGH HIGH",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ggt_ctcv6,
    compare = actual_ggt,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "BNRIND", "AVALU")
  )
})

### Haptoglobin decreased (NCICTCv4)
# Same as NCICTCv5
### Grade 1: <LLN

expected_hapt <- tibble::tribble(
  ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRL,
  "Not a term",            9,      10,     40,     NA_character_,  NA,
  NA_character_,           10,     10,     40,     NA_character_,  NA,
  "Haptoglobin decreased", 9,      10,     40,     NA_character_,  "1",
  "Haptoglobin decreased", 10,     10,     40,     NA_character_,  "0",
  "Haptoglobin decreased", 11,     10,     40,     NA_character_,  "0",
  # ANRHI missing - cannot grade
  "Haptoglobin decreased", NA,     10,     40,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Haptoglobin decreased", 10,     NA,     NA,     NA_character_,  NA,
)

input_hapt <- expected_hapt %>%
  select(-ATOXGRL)

## Test 31: CTCAEv4 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 31: CTCAEv4 Haptoglobin decreased", {
  actual_hapt <- derive_var_atoxgr_dir(
    input_hapt,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hapt,
    compare = actual_hapt,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 32a: CTCAEv5 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 32a: CTCAEv5 Haptoglobin decreased", {
  actual_hapt <- derive_var_atoxgr_dir(
    input_hapt,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hapt,
    compare = actual_hapt,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 32b: CTCAEv6 Haptoglobin decreased ----
test_that("derive_var_atoxgr_dir Test 32b: CTCAEv6 Haptoglobin decreased", {
  actual_hapt <- derive_var_atoxgr_dir(
    input_hapt,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hapt,
    compare = actual_hapt,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hemoglobin increased
# NCICTCAEv5 same as NCICTCAEv4 when BASE is normal
# SI unit is "g/L" USCV unit is "g/dL"
### Grade 3: Increase in >4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 2: Increase in >2 - 4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 1: Increase in >0 - 2 gm/dL above ULN or above baseline if baseline is above ULN

expected_hgbi_si <- tibble::tribble(
  ~ATOXDSCH,              ~AVAL, ~BASE, ~ANRHI, ~AVALU, ~ATOXGRH, ~TESTNUM, ~V5,
  "Not a term",           80,    120,   200,     "g/L",       NA,        1, "Y",
  NA_character_,          60,    50,    100,     "g/L",       NA,        2, "Y",
  # BASE greater than ANRHI
  "Hemoglobin increased", 106,   65,    60,      "g/L",      "3",        3, "N",
  "Hemoglobin increased", 105,   65,    60,      "g/L",      "2",        4, "N",
  "Hemoglobin increased", 86,    65,    60,      "g/L",      "2",        5, "N",
  "Hemoglobin increased", 85,    65,    60,      "g/L",      "1",        6, "N",
  "Hemoglobin increased", 66,    65,    60,      "g/L",      "1",        7, "N",
  "Hemoglobin increased", 65,    65,    60,      "g/L",      "0",        8, "N",
  "Hemoglobin increased", NA,    65,    60,      "g/L",       NA,        9, "N",
  # BASE less than or equal to ANRHI
  "Hemoglobin increased", 106,   60,    65,      "g/L",      "3",       10, "Y",
  "Hemoglobin increased", 105,   60,    65,      "g/L",      "2",       11, "Y",
  "Hemoglobin increased", 86,    60,    65,      "g/L",      "2",       12, "Y",
  "Hemoglobin increased", 85,    60,    65,      "g/L",      "1",       13, "Y",
  "Hemoglobin increased", 66,    60,    65,      "g/L",      "1",       14, "Y",
  "Hemoglobin increased", 65,    60,    65,      "g/L",      "0",       15, "Y",
  "Hemoglobin increased", NA,    60,    65,      "g/L",       NA,       16, "Y",
  # BASE missing
  "Hemoglobin increased", 106,   NA,    65,      "g/L",      "3",       17, "N",
  "Hemoglobin increased", 105,   NA,    65,      "g/L",      "2",       18, "N",
  "Hemoglobin increased", 86,    NA,    65,      "g/L",      "2",       19, "N",
  "Hemoglobin increased", 85,    NA,    65,      "g/L",      "1",       20, "N",
  "Hemoglobin increased", 66,    NA,    65,      "g/L",      "1",       21, "N",
  "Hemoglobin increased", 65,    NA,    65,      "g/L",      "0",       22, "N",
  "Hemoglobin increased", NA,    NA,    65,      "g/L",       NA,       23, "N",
  # Unit missing cannot grade
  "Hemoglobin increased", 200,   61,    65,         NA,       NA,       24, "Y",
  # ANRHI missing - cannot grade
  "Hemoglobin increased", 200,   60,    NA,      "g/L",       NA,       25, "Y",
  # AVAL missing cannot grade
  "Hemoglobin increased", NA,    60,    65,      "g/L",       NA,       26, "Y",
)

expected_hgbi_cv <- expected_hgbi_si %>%
  mutate(
    AVAL = AVAL / 10,
    BASE = BASE / 10,
    ANRHI = ANRHI / 10,
    AVALU = if_else(str_to_upper(AVALU) == "G/L", "g/dL", NA_character_)
  )

## Test 33a: CTCAEv4 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 33a: CTCAEv4 Hemoglobin increased (SI unit)", {
  input_hgbi_si <- expected_hgbi_si %>%
    select(-ATOXGRH)

  actual_hgbi_si <- derive_var_atoxgr_dir(
    input_hgbi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_si,
    compare = actual_hgbi_si,
    keys = c("TESTNUM")
  )
})

## Test 33b: CTCAEv4 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 33b: CTCAEv4 Hemoglobin increased (USCV unit)", {
  input_hgbi_cv <- expected_hgbi_cv %>%
    select(-ATOXGRH)

  actual_hgbi_cv <- derive_var_atoxgr_dir(
    input_hgbi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_cv,
    compare = actual_hgbi_cv,
    keys = c("TESTNUM")
  )
})

## Test 34a: CTCAEv5 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 34a: CTCAEv5 Hemoglobin increased (SI unit)", {
  expected_hgbi_si <- expected_hgbi_si %>%
    filter(V5 == "Y")

  input_hgbi_si <- expected_hgbi_si %>%
    select(-ATOXGRH)

  actual_hgbi_si <- derive_var_atoxgr_dir(
    input_hgbi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_si,
    compare = actual_hgbi_si,
    keys = c("TESTNUM")
  )
})

## Test 34b: CTCAEv5 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 34a: CTCAEv5 Hemoglobin increased (USCV unit)", {
  expected_hgbi_cv <- expected_hgbi_cv %>%
    filter(V5 == "Y")

  input_hgbi_cv <- expected_hgbi_cv %>%
    select(-ATOXGRH)

  actual_hgbi_cv <- derive_var_atoxgr_dir(
    input_hgbi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_cv,
    compare = actual_hgbi_cv,
    keys = c("TESTNUM")
  )
})

## Test 34c: CTCAEv6 Hemoglobin increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 34c: CTCAEv6 Hemoglobin increased (SI unit)", {
  expected_hgbi_si <- expected_hgbi_si %>%
    filter(V5 == "Y")

  input_hgbi_si <- expected_hgbi_si %>%
    select(-ATOXGRH)

  actual_hgbi_si <- derive_var_atoxgr_dir(
    input_hgbi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_si,
    compare = actual_hgbi_si,
    keys = c("TESTNUM")
  )
})

## Test 34d: CTCAEv6 Hemoglobin increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 34d: CTCAEv6 Hemoglobin increased (USCV unit)", {
  expected_hgbi_cv <- expected_hgbi_cv %>%
    filter(V5 == "Y")

  input_hgbi_cv <- expected_hgbi_cv %>%
    select(-ATOXGRH)

  actual_hgbi_cv <- derive_var_atoxgr_dir(
    input_hgbi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi_cv,
    compare = actual_hgbi_cv,
    keys = c("TESTNUM")
  )
})

### INR increased (NCICTCV4)
### NCICTCV5 different for grade 1
### Grade 3: >2.5 x ULN; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5 x ULN; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1 - 1.5 x ULN; >1 - 1.5 times above baseline if on anticoagulation

## Test 35: CTCAEv4 INR increased ----
test_that("derive_var_atoxgr_dir Test 35: CTCAEv4 INR increased", {
  expected_inri <- tibble::tribble(
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
  input_inri <- expected_inri %>%
    select(-ATOXGRH)

  actual_inri <- derive_var_atoxgr_dir(
    input_inri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_inri,
    compare = actual_inri,
    keys = c("TESTNUM")
  )
})

### INR increased (NCICTCV5)
### NCICTCV5 different to NCICTCv4 (do not use x ULN)
### Grade 3: >2.5; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1.2 - 1.5; >1 - 1.5 times above baseline if on anticoagulation

## Test 36a: CTCAEv5 INR increased ----
test_that("derive_var_atoxgr_dir Test 36a: CTCAEv5 INR increased", {
  expected_inri <- tibble::tribble(
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
  input_inri <- expected_inri %>%
    select(-ATOXGRH)

  actual_inri <- derive_var_atoxgr_dir(
    input_inri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_inri,
    compare = actual_inri,
    keys = c("TESTNUM")
  )
})

## Test 36b: CTCAEv6 INR increased ----
test_that("derive_var_atoxgr_dir Test 36b: CTCAEv6 INR increased", {
  expected_inri <- tibble::tribble(
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
  input_inri <- expected_inri %>%
    select(-ATOXGRH)

  actual_inri <- derive_var_atoxgr_dir(
    input_inri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_inri,
    compare = actual_inri,
    keys = c("TESTNUM")
  )
})

### Lipase increased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

expected_lip <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",       80,     120,    200,    NA_character_,  NA,
  NA_character_,      60,     50,     100,    NA_character_,  NA,
  "Lipase IncreaSed", 501,    0,      100,    NA_character_,  "4",
  "Lipase Increased", 500,    0,      100,    NA_character_,  "3",
  "Lipase Increased", 201,    0,      100,    NA_character_,  "3",
  "Lipase Increased", 200,    0,      100,    NA_character_,  "2",
  "Lipase Increased", 151,    0,      100,    NA_character_,  "2",
  "Lipase Increased", 150,    0,      100,    NA_character_,  "1",
  "Lipase Increased", 101,    0,      100,    NA_character_,  "1",
  "Lipase Increased", 100,    0,      100,    NA_character_,  "0",
  # ANRHI missing cannot grade
  "Lipase Increased", 200,    0,      NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Lipase Increased", NA,     0,      100,    NA_character_,  NA,
)
input_lip <- expected_lip %>%
  select(-ATOXGRH)

## Test 37: CTCAEv4 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 37: CTCAEv4 Lipase increased", {
  actual_lip <- derive_var_atoxgr_dir(
    input_lip,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lip,
    compare = actual_lip,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 38a: CTCAEv5 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 38a: CTCAEv5 Lipase increased", {
  actual_lip <- derive_var_atoxgr_dir(
    input_lip,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lip,
    compare = actual_lip,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Lipase increased NCICTCAEv6 (small difference to v4 and v5)
### Grade 4: >5.0 x ULN
### Grade 3: >3.0 - 5.0 x ULN
### Grade 2: >1.5 - 3.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

expected_lipv6 <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",       80,     120,    200,    NA_character_,  NA,
  NA_character_,      60,     50,     100,    NA_character_,  NA,
  "Lipase IncreaSed", 501,    0,      100,    NA_character_,  "4",
  "Lipase Increased", 500,    0,      100,    NA_character_,  "3",
  "Lipase Increased", 301,    0,      100,    NA_character_,  "3",
  "Lipase Increased", 300,    0,      100,    NA_character_,  "2",
  "Lipase Increased", 151,    0,      100,    NA_character_,  "2",
  "Lipase Increased", 150,    0,      100,    NA_character_,  "1",
  "Lipase Increased", 101,    0,      100,    NA_character_,  "1",
  "Lipase Increased", 100,    0,      100,    NA_character_,  "0",
  # ANRHI missing cannot grade
  "Lipase Increased", 200,    0,      NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Lipase Increased", NA,     0,      100,    NA_character_,  NA,
)
input_lipv6 <- expected_lipv6 %>%
  select(-ATOXGRH)


## Test 38b: CTCAEv6 Lipase increased ----
test_that("derive_var_atoxgr_dir Test 38b: CTCAEv6 Lipase increased", {
  actual_lipv6 <- derive_var_atoxgr_dir(
    input_lipv6,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lipv6,
    compare = actual_lipv6,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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

expected_lymd_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  # ANRLO not missing
  "Lymphocyte count decreased", 0.19,   140,    NA,     "10^9/L",  "4",      3,
  "Lymphocyte count decreased", 0.2,    140,    NA,     "10^9/L",  "3",      4,
  "Lymphocyte count decreased", 0.49,   140,    NA,     "10^9/L",  "3",      5,
  "Lymphocyte count decreased", 0.5,    140,    NA,     "10^9/L",  "2",      6,
  "Lymphocyte count decreased", 0.79,   140,    NA,     "10^9/L",  "2",      7,
  "Lymphocyte count decreased", 0.8,    140,    NA,     "10^9/L",  "1",      8,
  "Lymphocyte count decreased", 1.09,   1.1,    NA,     "10^9/L",  "1",      9,
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     "10^9/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Lymphocyte count decreased", 0.19,   NA,     NA,     "10^9/L",  "4",      11,
  "Lymphocyte count decreased", 0.2,    NA,     NA,     "10^9/L",  "3",      12,
  "Lymphocyte count decreased", 0.49,   NA,     NA,     "10^9/L",  "3",      13,
  "Lymphocyte count decreased", 0.5,    NA,     NA,     "10^9/L",  "2",      14,
  "Lymphocyte count decreased", 0.79,   NA,     NA,     "10^9/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Lymphocyte count decreased", 0.8,    NA,     NA,     "10^9/L",  NA,       16,
  "Lymphocyte count decreased", 1.09,   NA,     NA,     "10^9/L",  NA,       17,
  "Lymphocyte count decreased", 1.1,    NA,     NA,     "10^9/L",  NA,       18,
  # Unit missing cannot grade
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     "10^9/L",  "0",      20,
)
input_lymd_si <- expected_lymd_si %>%
  select(-ATOXGRL)

## Test 39a: CTCAEv4 Lymphocyte count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 39a: CTCAEv4 Lymphocyte count decreased (SI unit)", {
  actual_lymd_si <- derive_var_atoxgr_dir(
    input_lymd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_si,
    compare = actual_lymd_si,
    keys = c("TESTNUM")
  )
})

expected_lymd_cv <- expected_lymd_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )
input_lymd_cv <- expected_lymd_cv %>%
  select(-ATOXGRL)

## Test 39b: CTCAEv4 Lymphocyte count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 39a: CTCAEv4 Lymphocyte count decreased (USCV unit)", {
  actual_lymd_cv <- derive_var_atoxgr_dir(
    input_lymd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_cv,
    compare = actual_lymd_cv,
    keys = c("TESTNUM")
  )
})

expected_lymd_cv2 <- expected_lymd_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
input_lymd_cv2 <- expected_lymd_cv2 %>%
  select(-ATOXGRL)

## Test 39c: CTCAEv4 Lymphocyte count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 39c: CTCAEv4 Lymphocyte count decreased (legacy USCV unit)", {
  actual_lymd_cv2 <- derive_var_atoxgr_dir(
    input_lymd_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_cv2,
    compare = actual_lymd_cv2,
    keys = c("TESTNUM")
  )
})

## Test 40a: CTCAEv5 Lymphocyte count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 40a: CTCAEv5 Lymphocyte count decreased (SI unit)", {
  actual_lymd_si <- derive_var_atoxgr_dir(
    input_lymd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_si,
    compare = actual_lymd_si,
    keys = c("TESTNUM")
  )
})

## Test 40b: CTCAEv5 Lymphocyte count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 40b: CTCAEv5 Lymphocyte count decreased (USCV unit)", {
  actual_lymd_cv <- derive_var_atoxgr_dir(
    input_lymd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_cv,
    compare = actual_lymd_cv,
    keys = c("TESTNUM")
  )
})

expected_lymd_cv2 <- expected_lymd_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
input_lymd_cv2 <- expected_lymd_cv2 %>%
  select(-ATOXGRL)

## Test 40c: CTCAEv5 Lymphocyte count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 40c: CTCAEv5 Lymphocyte count decreased (legacy USCV unit)", {
  actual_lymd_cv2 <- derive_var_atoxgr_dir(
    input_lymd_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd_cv2,
    compare = actual_lymd_cv2,
    keys = c("TESTNUM")
  )
})

### Lymphocyte count increased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 3: >20,000/mm3
### Grade 2: >4000/mm3 - 20,000/mm3

expected_lymi_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  "Lymphocyte count increased", 21,     NA,     NA,     "10^9/L",  "3",      3,
  "Lymphocyte count increased", 20,     NA,     NA,     "10^9/L",  "2",      4,
  "Lymphocyte count increased", 4.1,    NA,     NA,     "10^9/L",  "2",      5,
  "Lymphocyte count increased", 4,      NA,     NA,     "10^9/L",  "0",      6,
  # Unit missing cannot grade
  "Lymphocyte count increased", 4,      NA,     NA,     NA,        NA,       7,
  # AVAL missing cannot grade
  "Lymphocyte count increased", NA,     NA,     NA,     "10^9/L",  NA,       8,
)
input_lymi_si <- expected_lymi_si %>%
  select(-ATOXGRH)

## Test 41a: CTCAEv4 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 41a: CTCAEv4 Lymphocyte count increased (SI unit)", {
  actual_lymi_si <- derive_var_atoxgr_dir(
    input_lymi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_si,
    compare = actual_lymi_si,
    keys = c("TESTNUM")
  )
})

expected_lymi_cv1 <- expected_lymi_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

expected_lymi_cv2 <- expected_lymi_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 8
  )

expected_lymi_cv <- bind_rows(expected_lymi_cv1, expected_lymi_cv2)

input_lymi_cv <- expected_lymi_cv %>%
  select(-ATOXGRH)

## Test 41b: CTCAEv4 Lymphocyte count increased  (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 41b: CTCAEv4 Lymphocyte count increased  (USCV unit)", {
  actual_lymi_cv <- derive_var_atoxgr_dir(
    input_lymi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_cv,
    compare = actual_lymi_cv,
    keys = c("TESTNUM")
  )
})


## Test 42a: CTCAEv5 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 42a: CTCAEv5 Lymphocyte count increased (SI unit)", {
  actual_lymi_si <- derive_var_atoxgr_dir(
    input_lymi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_si,
    compare = actual_lymi_si,
    keys = c("TESTNUM")
  )
})

## Test 42b: CTCAEv5 Lymphocyte count increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 42b: CTCAEv5 Lymphocyte count increased (USCV unit)", {
  actual_lymi_cv <- derive_var_atoxgr_dir(
    input_lymi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_cv,
    compare = actual_lymi_cv,
    keys = c("TESTNUM")
  )
})

## Test 42c: CTCAEv6 Lymphocyte count increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 42c: CTCAEv6 Lymphocyte count increased (SI unit)", {
  actual_lymi_si <- derive_var_atoxgr_dir(
    input_lymi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_si,
    compare = actual_lymi_si,
    keys = c("TESTNUM")
  )
})

## Test 42d: CTCAEv6 Lymphocyte count increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 42d: CTCAEv6 Lymphocyte count increased (USCV unit)", {
  actual_lymi_cv <- derive_var_atoxgr_dir(
    input_lymi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi_cv,
    compare = actual_lymi_cv,
    keys = c("TESTNUM")
  )
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

expected_neut_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  # ANRLO not missing
  "Neutrophil count decreased", 0.49,   2,      NA,     "10^9/L",  "4",      3,
  "Neutrophil count decreased", 0.5,    2,      NA,     "10^9/L",  "3",      4,
  "Neutrophil count decreased", 0.99,   2,      NA,     "10^9/L",  "3",      5,
  "Neutrophil count decreased", 1,      2,      NA,     "10^9/L",  "2",      6,
  "Neutrophil count decreased", 1.49,   2,      NA,     "10^9/L",  "2",      7,
  "Neutrophil count decreased", 1.5,    2,      NA,     "10^9/L",  "1",      8,
  "Neutrophil count decreased", 1.9,    2,      NA,     "10^9/L",  "1",      9,
  "Neutrophil count decreased", 2,      2,      NA,     "10^9/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Neutrophil count decreased", 0.49,   NA,     NA,     "10^9/L",  "4",      11,
  "Neutrophil count decreased", 0.5,    NA,     NA,     "10^9/L",  "3",      12,
  "Neutrophil count decreased", 0.99,   NA,     NA,     "10^9/L",  "3",      13,
  "Neutrophil count decreased", 1,      NA,     NA,     "10^9/L",  "2",      14,
  "Neutrophil count decreased", 1.49,   NA,     NA,     "10^9/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Neutrophil count decreased", 1.5,    NA,     NA,     "10^9/L",  NA,       16,
  "Neutrophil count decreased", 1.9,    NA,     NA,     "10^9/L",  NA,       17,
  "Neutrophil count decreased", 2,      NA,     NA,     "10^9/L",  NA,       18,
  # Unit missing cannot grade
  "Neutrophil count decreased", 2,      2,      NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Neutrophil count decreased", NA,     2,      NA,     "10^9/L",  NA,       20,
)
input_neut_si <- expected_neut_si %>%
  select(-ATOXGRL)


## Test 43a: CTCAEv4 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 43a: CTCAEv4 Neutrophil count decreased (SI unit)", {
  actual_neut_si <- derive_var_atoxgr_dir(
    input_neut_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_si,
    compare = actual_neut_si,
    keys = c("TESTNUM")
  )
})

expected_neut_cv1 <- expected_neut_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

expected_neut_cv2 <- expected_neut_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 20
  )

expected_neut_cv <- bind_rows(expected_neut_cv1, expected_neut_cv2)

input_neut_cv <- expected_neut_cv %>%
  select(-ATOXGRL)

## Test 43b: CTCAEv4 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 43b: CTCAEv4 Neutrophil count decreased (USCV unit)", {
  actual_neut_cv <- derive_var_atoxgr_dir(
    input_neut_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_cv,
    compare = actual_neut_cv,
    keys = c("TESTNUM")
  )
})


## Test 44a: CTCAEv5 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 44a: CTCAEv5 Neutrophil count decreased (SI unit)", {
  actual_neut_si <- derive_var_atoxgr_dir(
    input_neut_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_si,
    compare = actual_neut_si,
    keys = c("TESTNUM")
  )
})

## Test 44b: CTCAEv5 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44b: CTCAEv5 Neutrophil count decreased (USCV unit)", {
  actual_neut_cv <- derive_var_atoxgr_dir(
    input_neut_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_cv,
    compare = actual_neut_cv,
    keys = c("TESTNUM")
  )
})

### Neutrophil count decreased NCICTCAEv6 (different to v4 and v5)
### SI unit is 10^9/L
### CV unit is 10^3/uL (= SI unit)
### Legacy CV unit is 10^3/mL (= 1000 * SI unit)
### Grade 4: <0.1 x 10^9/L
### Grade 3: <0.5 - 0.1 x 10^9/L
### Grade 2: <1.0 - 0.5 x 10^9/L
### Grade 1: <1.5 - 1.0 x 10^9/L

expected_neut_siv6 <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     "10^9/L",  NA,       1,
  NA_character_,                60,     "10^9/L",  NA,       2,
  # ANRLO not missing
  "Neutrophil count decreased", 0.09,   "10^9/L",  "4",      3,
  "Neutrophil count decreased", 0.1,    "10^9/L",  "3",      4,
  "Neutrophil count decreased", 0.49,   "10^9/L",  "3",      5,
  "Neutrophil count decreased", 0.5,    "10^9/L",  "2",      6,
  "Neutrophil count decreased", 0.99,   "10^9/L",  "2",      7,
  "Neutrophil count decreased", 1,      "10^9/L",  "1",      8,
  "Neutrophil count decreased", 1.49,   "10^9/L",  "1",      9,
  "Neutrophil count decreased", 1.5,    "10^9/L",  "0",      10,
  # Unit missing cannot grade
  "Neutrophil count decreased", NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Neutrophil count decreased", NA,     "10^9/L",  NA,       20,
)
input_neut_siv6 <- expected_neut_siv6 %>%
  select(-ATOXGRL)

## Test 44c: CTCAEv6 Neutrophil count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 44c: CTCAEv6 Neutrophil count decreased (SI unit)", {
  actual_neut_siv6 <- derive_var_atoxgr_dir(
    input_neut_siv6,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_siv6,
    compare = actual_neut_siv6,
    keys = c("TESTNUM")
  )
})

## Test 44d: CTCAEv6 Neutrophil count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44d: CTCAEv6 Neutrophil count decreased (USCV unit)", {
  expected_neut_uscv_v6 <- expected_neut_siv6 %>%
    mutate(
      AVALU = if_else(AVALU == "10^9/L", "10^3/uL", AVALU)
    )

  input_neut_uscv_v6 <- expected_neut_uscv_v6 %>%
    select(-ATOXGRL)

  actual_neut_uscv_v6 <- derive_var_atoxgr_dir(
    input_neut_uscv_v6,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_uscv_v6,
    compare = actual_neut_uscv_v6,
    keys = c("TESTNUM")
  )
})

## Test 44e: CTCAEv6 Neutrophil count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 44e: CTCAEv6 Neutrophil count decreased (legacy USCV unit)", {
  expected_neut_uscv2_v6 <- expected_neut_siv6 %>%
    mutate(
      AVALU = if_else(AVALU == "10^9/L", "10^3/mL", AVALU),
      AVAL = AVAL * 1000
    )

  input_neut_uscv2_v6 <- expected_neut_uscv2_v6 %>%
    select(-ATOXGRL)

  actual_neut_uscv2_v6 <- derive_var_atoxgr_dir(
    input_neut_uscv2_v6,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut_uscv2_v6,
    compare = actual_neut_uscv2_v6,
    keys = c("TESTNUM")
  )
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

expected_plate_si <- tibble::tribble(
  ~ATOXDSCL,                  ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",               80,    120,    200,    "10^9/L", NA,       1,
  NA_character_,              60,    50,     100,    "10^9/L", NA,       2,
  # ANRLO not missing
  "Platelet count decreased", 24,    100,    150,    "10^9/L", "4",      3,
  "Platelet count decreased", 25,    100,    150,    "10^9/L", "3",      4,
  "Platelet count decreased", 49,    100,    150,    "10^9/L", "3",      5,
  "Platelet count decreased", 50,    100,    150,    "10^9/L", "2",      6,
  "Platelet count decreased", 74,    100,    150,    "10^9/L", "2",      7,
  "Platelet count decreased", 75,    100,    150,    "10^9/L", "1",      8,
  "Platelet count decreased", 99,    100,    NA,     "10^9/L", "1",      9,
  "Platelet count decreased", 100,   100,    NA,     "10^9/L", "0",      10,
  # ANRLO missing - can grade 2-4      1,
  "Platelet count decreased", 24,    NA,     NA,     "10^9/L", "4",      11,
  "Platelet count decreased", 25,    NA,     NA,     "10^9/L", "3",      12,
  "Platelet count decreased", 49,    NA,     NA,     "10^9/L", "3",      13,
  "Platelet count decreased", 50,    NA,     NA,     "10^9/L", "2",      14,
  "Platelet count decreased", 74,    NA,     NA,     "10^9/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Platelet count decreased", 75,    NA,     NA,     "10^9/L", NA,       16,
  "Platelet count decreased", 99,    NA,     NA,     "10^9/L", NA,       17,
  "Platelet count decreased", 100,   NA,     NA,     "10^9/L", NA,       18,
  # Unit missing cannot grade
  "Platelet count decreased", 100,   100,    NA,     NA,       NA,       19,
  # AVAL missing cannot grade
  "Platelet count decreased", NA,    100,    NA,     "10^9/L", NA,       20,
)
input_plate_si <- expected_plate_si %>%
  select(-ATOXGRL)

## Test 45a: CTCAEv4 Platelet count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 45a: CTCAEv4 Platelet count decreased (SI unit)", {
  actual_plate_si <- derive_var_atoxgr_dir(
    input_plate_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_si,
    compare = actual_plate_si,
    keys = c("TESTNUM")
  )
})

expected_plate_cv <- expected_plate_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )
input_plate_cv <- expected_plate_cv %>%
  select(-ATOXGRL)

## Test 45b: CTCAEv4 Platelet count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 45b: CTCAEv4 Platelet count decreased (USCV unit)", {
  actual_plate_cv <- derive_var_atoxgr_dir(
    input_plate_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_cv,
    compare = actual_plate_cv,
    keys = c("TESTNUM")
  )
})

expected_plate_cv2 <- expected_plate_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
input_plate_cv2 <- expected_plate_cv2 %>%
  select(-ATOXGRL)

## Test 45c: CTCAEv4 Platelet count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 45c: CTCAEv4 Platelet count decreased (legacy USCV unit)", {
  actual_plate_cv2 <- derive_var_atoxgr_dir(
    input_plate_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_cv2,
    compare = actual_plate_cv2,
    keys = c("TESTNUM")
  )
})

## Test 46a: CTCAEv5 Platelet count decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 46a: CTCAEv5 Platelet count decreased (SI unit)", {
  actual_plate_si <- derive_var_atoxgr_dir(
    input_plate_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_si,
    compare = actual_plate_si,
    keys = c("TESTNUM")
  )
})

## Test 46b: CTCAEv5 Platelet count decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46b: CTCAEv5 Platelet count decreased (USCV unit)", {
  actual_plate_cv <- derive_var_atoxgr_dir(
    input_plate_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_cv,
    compare = actual_plate_cv,
    keys = c("TESTNUM")
  )
})

## Test 46c: CTCAEv5 Platelet count decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46c: CTCAEv5 Platelet count decreased (legacy USCV unit)", {
  actual_plate_cv2 <- derive_var_atoxgr_dir(
    input_plate_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_cv2,
    compare = actual_plate_cv2,
    keys = c("TESTNUM")
  )
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

expected_plate_v6_si <- expected_plate_si %>%
  mutate(
    ATOXDSCL = if_else(
      ATOXDSCL == "Platelet count decreased", "Thrombocytopenia", ATOXDSCL
    ),
    AVAL = case_when(
      AVAL == 24 ~ 9,
      AVAL == 25 ~ 10,
      TRUE ~ AVAL
    )
  )

## Test 46d: CTCAEv6 Thrombocytopenia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 46d: CTCAEv6 Thrombocytopenia (SI unit)", {
  input_plate_v6_si <- expected_plate_v6_si %>%
    select(-ATOXGRL)

  actual_plate_si <- derive_var_atoxgr_dir(
    input_plate_v6_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_v6_si,
    compare = actual_plate_si,
    keys = c("TESTNUM")
  )
})

## Test 46e: CTCAEv6 Thrombocytopenia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46e: CTCAEv6 Thrombocytopenia (USCV unit)", {
  expected_plate_v6_cv <- expected_plate_v6_si %>%
    mutate(AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU))

  input_plate_v6_cv <- expected_plate_v6_cv %>%
    select(-ATOXGRL)

  actual_plate_cv <- derive_var_atoxgr_dir(
    input_plate_v6_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_v6_cv,
    compare = actual_plate_cv,
    keys = c("TESTNUM")
  )
})

## Test 46f: CTCAEv6 Thrombocytopenia (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 46e: CTCAEv6 Platelet count decreased (legacy USCV unit)", {
  expected_plate_v6_cv2 <- expected_plate_v6_si %>%
    mutate(
      AVAL = AVAL * 1000,
      ANRLO = ANRLO * 1000,
      ANRHI = ANRHI * 1000,
      AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
    )

  input_plate_v6_cv2 <- expected_plate_v6_cv2 %>%
    select(-ATOXGRL)

  actual_plate_cv2 <- derive_var_atoxgr_dir(
    input_plate_v6_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate_v6_cv2,
    compare = actual_plate_cv2,
    keys = c("TESTNUM")
  )
})

### Serum amylase increased
### NCICTCAEv4 and NCICTCAEv5 criteria essentially the same
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

expected_seri <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRH,
  "Not a term",              80,     120,    200,    NA_character_,  NA,
  NA_character_,             60,     50,     100,    NA_character_,  NA,
  "Serum amylase increased", 501,    0,      100,    NA_character_,  "4",
  "Serum amylase increased", 500,    0,      100,    NA_character_,  "3",
  "Serum amylase increased", 201,    0,      100,    NA_character_,  "3",
  "Serum amylase increased", 200,    0,      100,    NA_character_,  "2",
  "Serum amylase increased", 151,    0,      100,    NA_character_,  "2",
  "Serum amylase increased", 150,    0,      100,    NA_character_,  "1",
  "Serum amylase increased", 101,    0,      100,    NA_character_,  "1",
  "Serum amylase increased", 100,    0,      100,    NA_character_,  "0",
  # ANRHI missing cannot grade
  "Serum amylase increased", 200,    0,      NA,     NA_character_,  NA,
  # AVAL missing cannot grade
  "Serum amylase increased", NA,     0,      100,    NA_character_,  NA,
)
input_seri <- expected_seri %>%
  select(-ATOXGRH)

## Test 47a: CTCAEv4 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 47a: CTCAEv4 Serum amylase increased (SI unit)", {
  actual_seri <- derive_var_atoxgr_dir(
    input_seri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_seri,
    compare = actual_seri,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 47b: CTCAEv4 Serum amylase increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 47b: CTCAEv4 Serum amylase increased (USCV unit)", {
  actual_seri <- derive_var_atoxgr_dir(
    input_seri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_seri,
    compare = actual_seri,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 48a: CTCAEv5 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 48a: CTCAEv5 Serum amylase increased (SI unit)", {
  actual_seri <- derive_var_atoxgr_dir(
    input_seri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_seri,
    compare = actual_seri,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 48b: CTCAEv5 Serum amylase increased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 48b: CTCAEv5 Serum amylase increased (USCV unit)", {
  actual_seri <- derive_var_atoxgr_dir(
    input_seri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_seri,
    compare = actual_seri,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 48c: CTCAEv6 Serum amylase increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 48c: CTCAEv6 Serum amylase increased (SI unit)", {
  actual_seri <- derive_var_atoxgr_dir(
    input_seri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_seri,
    compare = actual_seri,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### White blood cell decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <1.0 x 10e9 /L
### Grade 3: <2.0 - 1.0 x 10e9 /L
### Grade 2: <3.0 - 2.0 x 10e9 /L
### Grade 1: <LLN - 3.0 x 10e9 /L

expected_wbcd <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",                 1,     5,      15,     "10^9/L", NA,       1,
  NA_character_,                2,     5,      15,     "10^9/L", NA,       2,
  # ANRLO not missing
  "White blood cell decreased", 0.9,   5,      15,     "10^9/L", "4",      3,
  "White blood cell decreased", 1,     5,      15,     "10^9/L", "3",      4,
  "White blood cell decreased", 1.9,   5,      15,     "10^9/L", "3",      5,
  "White blood cell decreased", 2,     5,      15,     "10^9/L", "2",      6,
  "White blood cell decreased", 2.9,   5,      15,     "10^9/L", "2",      7,
  "White blood cell decreased", 3,     5,      15,     "10^9/L", "1",      8,
  "White blood cell decreased", 4.9,   5,      15,     "10^9/L", "1",      9,
  "White blood cell decreased", 5,     5,      15,     "10^9/L", "0",      10,
  # ANRLO missing - can grade 2-4
  "White blood cell decreased", 0.9,   NA,     NA,     "10^9/L", "4",      11,
  "White blood cell decreased", 1,     NA,     NA,     "10^9/L", "3",      12,
  "White blood cell decreased", 1.9,   NA,     NA,     "10^9/L", "3",      13,
  "White blood cell decreased", 2,     NA,     NA,     "10^9/L", "2",      14,
  "White blood cell decreased", 2.9,   NA,     NA,     "10^9/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       16,
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       17,
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       18,
  # Unit missing cannot grade
  "White blood cell decreased", 3,     100,    NA,     NA,       NA,       19,
  # AVAL missing cannot grade
  "White blood cell decreased", NA,    100,    NA,     "10^9/L", NA,       20,
)
input_wbcd <- expected_wbcd %>%
  select(-ATOXGRL)

expected_wbcd_uscv1 <- expected_wbcd %>%
  mutate(AVALU = if_else(AVALU == "10^9/L", "10^3/uL", AVALU))

expected_wbcd_uscv2 <- expected_wbcd %>%
  mutate(
    ANRLO = if_else(!is.na(ANRLO), ANRLO * 1000, NA),
    AVAL = if_else(!is.na(AVAL), AVAL * 1000, NA),
    AVALU = if_else(AVALU == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 20
  )

expected_wbcd_uscv <- bind_rows(expected_wbcd_uscv1, expected_wbcd_uscv2)

input_wbcd_uscv <- expected_wbcd_uscv %>%
  select(-ATOXGRL)


## Test 49a: CTCAEv4 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 49a: CTCAEv4 White blood cell decreased (SI unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
})

## Test 49b: CTCAEv4 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 49b: CTCAEv4 White blood cell decreased (USCV unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd_uscv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_uscv,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
})


## Test 50a: CTCAEv5 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 50a: CTCAEv5 White blood cell decreased (SI unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
})

## Test 50b: CTCAEv5 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 50b: CTCAEv5 White blood cell decreased (USCV unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd_uscv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_uscv,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
})

## Test 50c: CTCAEv6 White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 50c: CTCAEv6 White blood cell decreased (SI unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
})

## Test 50d: CTCAEv6 White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 50d: CTCAEv6 White blood cell decreased (USCV unit)", {
  actual_wbcd <- derive_var_atoxgr_dir(
    input_wbcd_uscv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_uscv,
    compare = actual_wbcd,
    keys = c("TESTNUM")
  )
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

expected_calci_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    0,      2.5,    "mmol/L",  NA,       1,
  NA_character_,   3.5,    0,      2.5,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypercalcemia", 3.5,    0,      2.5,    "mmol/L",  "4",      3,
  "Hypercalcemia", 3.4,    0,      2.5,    "mmol/L",  "3",      4,
  "Hypercalcemia", 3.2,    0,      2.5,    "mmol/L",  "3",      5,
  "Hypercalcemia", 3.1,    0,      2.5,    "mmol/L",  "2",      6,
  "Hypercalcemia", 3,      0,      2.5,    "mmol/L",  "2",      7,
  "Hypercalcemia", 2.9,    0,      2.5,    "mmol/L",  "1",      8,
  "Hypercalcemia", 2.6,    0,      2.5,    "mmol/L",  "1",      9,
  "Hypercalcemia", 2.5,    0,      2.5,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia", 3.5,    0,      NA,     "mmol/L",  "4",      11,
  "Hypercalcemia", 3.4,    0,      NA,     "mmol/L",  "3",      12,
  "Hypercalcemia", 3.2,    0,      NA,     "mmol/L",  "3",      13,
  "Hypercalcemia", 3.1,    0,      NA,     "mmol/L",  "2",      14,
  "Hypercalcemia", 3,      0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypercalcemia", 2.9,    0,      NA,     "mmol/L",  NA,       16,
  "Hypercalcemia", 2.6,    0,      NA,     "mmol/L",  NA,       17,
  "Hypercalcemia", 2.5,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypercalcemia", 2.5,    0,      2.5,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia", NA,     0,      2.5,    "mmol/L",  NA,       20,
)
input_calci_si <- expected_calci_si %>%
  select(-ATOXGRH)

expected_calci_cv <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI,   ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",    13.6,   0,      10.5,    "mg/dL",         NA,        1,
  NA_character_,   13.6,   0,      10.5,    "mg/dL",         NA,        2,
  # ANRHI not missing
  "Hypercalcemia", 13.6,   0,      10.5,    "mg/dL",        "4",        3,
  "Hypercalcemia", 13.5,   0,      10.5,    "mg/dL",        "3",        4,
  "Hypercalcemia", 12.6,   0,      10.5,    "mg/dL",        "3",        5,
  "Hypercalcemia", 12.5,   0,      10.5,    "mg/dL",        "2",        6,
  "Hypercalcemia", 11.6,   0,      10.5,    "mg/dL",        "2",        7,
  "Hypercalcemia", 11.5,   0,      10.5,    "mg/dL",        "1",        8,
  "Hypercalcemia", 10.6,   0,      10.5,    "mg/dL",        "1",        9,
  "Hypercalcemia", 10.5,   0,      10.5,    "mg/dL",        "0",       10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia", 13.6,   0,      NA,      "mg/dL",        "4",       11,
  "Hypercalcemia", 13.5,   0,      NA,      "mg/dL",        "3",       12,
  "Hypercalcemia", 12.6,   0,      NA,      "mg/dL",        "3",       13,
  "Hypercalcemia", 12.5,   0,      NA,      "mg/dL",        "2",       14,
  "Hypercalcemia", 11.6,   0,      NA,      "mg/dL",        "2",       15,
  # ANRHI missing - can NOT grade 0 or 1mg/d
  "Hypercalcemia", 11.5,   0,      NA,      "mg/dL",         NA,       16,
  "Hypercalcemia", 10.6,   0,      NA,      "mg/dL",         NA,       17,
  "Hypercalcemia", 10.5,   0,      NA,      "mg/dL",         NA,       18,
  # Unit missing cannot grade
  "Hypercalcemia", 10.5,   0,      10.5,         NA,         NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia", NA,     0,      10.5,    "mg/dL",         NA,       20,
)
input_calci_cv <- expected_calci_cv %>%
  select(-ATOXGRH)

## Test 51a: CTCAEv4 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 51a: CTCAEv4 Hypercalcemia (SI unit)", {
  actual_calci_si <- derive_var_atoxgr_dir(
    input_calci_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_si,
    compare = actual_calci_si,
    keys = c("TESTNUM")
  )
})

## Test 51b: CTCAEv4 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 51b: CTCAEv4 Hypercalcemia (USCV unit)", {
  actual_calci_cv <- derive_var_atoxgr_dir(
    input_calci_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_cv,
    compare = actual_calci_cv,
    keys = c("TESTNUM")
  )
})

## Test 52a: CTCAEv5 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 52a: CTCAEv5 Hypercalcemia (SI unit)", {
  actual_calci_si <- derive_var_atoxgr_dir(
    input_calci_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_si,
    compare = actual_calci_si,
    keys = c("TESTNUM")
  )
})

## Test 52b: CTCAEv5 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 52b: CTCAEv5 Hypercalcemia (USCV unit)", {
  actual_calci_cv <- derive_var_atoxgr_dir(
    input_calci_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_cv,
    compare = actual_calci_cv,
    keys = c("TESTNUM")
  )
})

## Test 52c: CTCAEv6 Hypercalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 52c: CTCAEv6 Hypercalcemia (SI unit)", {
  actual_calci_si <- derive_var_atoxgr_dir(
    input_calci_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_si,
    compare = actual_calci_si,
    keys = c("TESTNUM")
  )
})

## Test 52c: CTCAEv6 Hypercalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 52c: CTCAEv6 Hypercalcemia (USCV unit)", {
  actual_calci_cv <- derive_var_atoxgr_dir(
    input_calci_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_cv,
    compare = actual_calci_cv,
    keys = c("TESTNUM")
  )
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

expected_calioni_si <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",              1.9,    0,      1.3,    "mmol/L",  NA,       1,
  NA_character_,             1.9,    0,      1.3,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypercalcemia (Ionized)", 1.9,    0,      1.3,    "mmol/L",  "4",      3,
  "Hypercalcemia (Ionized)", 1.8,    0,      1.3,    "mmol/L",  "3",      4,
  "Hypercalcemia (Ionized)", 1.7,    0,      1.3,    "mmol/L",  "3",      5,
  "Hypercalcemia (Ionized)", 1.6,    0,      1.3,    "mmol/L",  "2",      6,
  "Hypercalcemia (Ionized)", 1.51,   0,      1.3,    "mmol/L",  "2",      7,
  "Hypercalcemia (Ionized)", 1.5,    0,      1.3,    "mmol/L",  "1",      8,
  "Hypercalcemia (Ionized)", 1.4,    0,      1.3,    "mmol/L",  "1",      9,
  "Hypercalcemia (Ionized)", 1.3,    0,      1.3,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia (Ionized)", 1.9,    0,      NA,     "mmol/L",  "4",      11,
  "Hypercalcemia (Ionized)", 1.8,    0,      NA,     "mmol/L",  "3",      12,
  "Hypercalcemia (Ionized)", 1.7,    0,      NA,     "mmol/L",  "3",      13,
  "Hypercalcemia (Ionized)", 1.6,    0,      NA,     "mmol/L",  "2",      14,
  "Hypercalcemia (Ionized)", 1.51,   0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypercalcemia (Ionized)", 1.5,    0,      NA,     "mmol/L",  NA,       16,
  "Hypercalcemia (Ionized)", 1.4,    0,      NA,     "mmol/L",  NA,       17,
  "Hypercalcemia (Ionized)", 1.3,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade       1,
  "Hypercalcemia (Ionized)", 1.3,    0,      1.3,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia (Ionized)", NA,     0,      1.3,    "mmol/L",  NA,       20,
)
input_calioni_si <- expected_calioni_si %>%
  select(-ATOXGRH)

expected_calioni_cv <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",               7.3,    0,       5,       "mg/dL",        NA,  1,
  NA_character_,              7.3,    0,       5,       "mg/dL",        NA,  2,
  # ANRHI not missing
  "Hypercalcemia (Ionized)",  7.3,    0,       5,       "mg/dL",       "4",  3,
  "Hypercalcemia (Ionized)",  7.2,    0,       5,       "mg/dL",       "3",  4,
  "Hypercalcemia (Ionized)",  6.5,    0,       5,       "mg/dL",       "3",  5,
  "Hypercalcemia (Ionized)",  6.4,    0,       5,       "mg/dL",       "2",  6,
  "Hypercalcemia (Ionized)",  6.1,    0,       5,       "mg/dL",       "2",  7,
  "Hypercalcemia (Ionized)",  6.0,    0,       5,       "mg/dL",       "1",  8,
  "Hypercalcemia (Ionized)",  5.1,    0,       5,       "mg/dL",       "1",  9,
  "Hypercalcemia (Ionized)",  5,      0,       5,       "mg/dL",       "0",  10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia (Ionized)",  7.3,    0,       NA,      "mg/dL",       "4",  11,
  "Hypercalcemia (Ionized)",  7.2,    0,       NA,      "mg/dL",       "3",  12,
  "Hypercalcemia (Ionized)",  6.5,    0,       NA,      "mg/dL",       "3",  13,
  "Hypercalcemia (Ionized)",  6.4,    0,       NA,      "mg/dL",       "2",  14,
  "Hypercalcemia (Ionized)",  6.1,    0,       NA,      "mg/dL",       "2",  15,
  # ANRHI missing - can NOT grade 0 or 1mg/d
  "Hypercalcemia (Ionized)",  6.0,    0,       NA,      "mg/dL",        NA,  16,
  "Hypercalcemia (Ionized)",  5.1,    0,       NA,      "mg/dL",        NA,  17,
  "Hypercalcemia (Ionized)",  5,      0,       NA,      "mg/dL",        NA,  18,
  # Unit missing cannot grade
  "Hypercalcemia (Ionized)",  5,      0,       5,            NA,        NA,  19,
  # AVAL missing cannot grade
  "Hypercalcemia (Ionized)",  NA,     0,       5,       "mg/dL",        NA,  20,
)
input_calioni_cv <- expected_calioni_cv %>%
  select(-ATOXGRH)

## Test 53a: CTCAEv4 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 53a: CTCAEv4 Hypercalcemia (Ionized) (SI unit)", {
  actual_calioni_si <- derive_var_atoxgr_dir(
    input_calioni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_si,
    compare = actual_calioni_si,
    keys = c("TESTNUM")
  )
})

## Test 53b: CTCAEv4 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 53b: CTCAEv4 Hypercalcemia (Ionized) (USCV unit)", {
  actual_calioni_cv <- derive_var_atoxgr_dir(
    input_calioni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_cv,
    compare = actual_calioni_cv,
    keys = c("TESTNUM")
  )
})

## Test 54a: CTCAEv5 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 54a: CTCAEv5 Hypercalcemia (Ionized) (SI unit)", {
  actual_calioni_si <- derive_var_atoxgr_dir(
    input_calioni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_si,
    compare = actual_calioni_si,
    keys = c("TESTNUM")
  )
})

## Test 54b: CTCAEv5 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 54b: CTCAEv5 Hypercalcemia (Ionized) (USCV unit)", {
  actual_calioni_cv <- derive_var_atoxgr_dir(
    input_calioni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_cv,
    compare = actual_calioni_cv,
    keys = c("TESTNUM")
  )
})

## Test 54c: CTCAEv6 Hypercalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 54c: CTCAEv6 Hypercalcemia (Ionized) (SI unit)", {
  actual_calioni_si <- derive_var_atoxgr_dir(
    input_calioni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_si,
    compare = actual_calioni_si,
    keys = c("TESTNUM")
  )
})

## Test 54d: CTCAEv6 Hypercalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 54d: CTCAEv6 Hypercalcemia (Ionized) (USCV unit)", {
  actual_calioni_cv <- derive_var_atoxgr_dir(
    input_calioni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_cv,
    compare = actual_calioni_cv,
    keys = c("TESTNUM")
  )
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

expected_glycfi <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",              27.9,   0,      5.3,    "mmol/L",  NA,       1,
  NA_character_,             27.9,   0,      5.3,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hyperglycemia (Fasting)", 27.9,   0,      5.3,    "mmol/L",  "4",      3,
  "Hyperglycemia (Fasting)", 27.8,   0,      5.3,    "mmol/L",  "3",      4,
  "Hyperglycemia (Fasting)", 14,     0,      5.3,    "mmol/L",  "3",      5,
  "Hyperglycemia (Fasting)", 13.9,   0,      5.3,    "mmol/L",  "2",      6,
  "Hyperglycemia (Fasting)", 9,      0,      5.3,    "mmol/L",  "2",      7,
  "Hyperglycemia (Fasting)", 8.9,    0,      5.3,    "mmol/L",  "1",      8,
  "Hyperglycemia (Fasting)", 5.4,    0,      5.3,    "mmol/L",  "1",      9,
  "Hyperglycemia (Fasting)", 5.3,    0,      5.3,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hyperglycemia (Fasting)", 27.9,   0,      NA,     "mmol/L",  "4",      11,
  "Hyperglycemia (Fasting)", 27.8,   0,      NA,     "mmol/L",  "3",      12,
  "Hyperglycemia (Fasting)", 14,     0,      NA,     "mmol/L",  "3",      13,
  "Hyperglycemia (Fasting)", 13.9,   0,      NA,     "mmol/L",  "2",      14,
  "Hyperglycemia (Fasting)", 9,      0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hyperglycemia (Fasting)", 8.9,    0,      NA,     "mmol/L",  NA,       16,
  "Hyperglycemia (Fasting)", 5.4,    0,      NA,     "mmol/L",  NA,       17,
  "Hyperglycemia (Fasting)", 5.3,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hyperglycemia (Fasting)", 5.3,    0,      5.3,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hyperglycemia (Fasting)", NA,     0,      5.3,    "mmol/L",  NA,       20,
)

input_glycfi <- expected_glycfi %>%
  select(-ATOXGRH)

## Test 55a: CTCAEv4 Hyperglycemia (Fasting) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 55a: CTCAEv4 Hyperglycemia (Fasting) (SI unit)", {
  actual_glycfi <- derive_var_atoxgr_dir(
    input_glycfi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycfi,
    compare = actual_glycfi,
    keys = c("TESTNUM")
  )
})

## Test 55b: CTCAEv4 Hyperglycemia (Fasting) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 55b: CTCAEv4 Hyperglycemia (Fasting) (USCV unit)", {
  expected_glycfi <- tibble::tribble(
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
  input_glycfi <- expected_glycfi %>%
    select(-ATOXGRH)

  actual_glycfi <- derive_var_atoxgr_dir(
    input_glycfi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycfi,
    compare = actual_glycfi,
    keys = c("TESTNUM")
  )
})

## Test 55c: CTCAEv6 Hyperglycemia (Fasting) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 55c: CTCAEv6 Hyperglycemia (Fasting) (SI unit)", {
  actual_glycfi <- derive_var_atoxgr_dir(
    input_glycfi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycfi,
    compare = actual_glycfi,
    keys = c("TESTNUM")
  )
})

### Hyperglycemia (NCICTCAEv4 and NCICTCAEv6)
### not included in NCICTCAEv5
### SI unit
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L

### CV unit
### Grade 4: >500 mg/dL
### Grade 3: >250 - 500 mg/dL

expected_glyci <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    27.9,   0,      5.3,    "mmol/L",  NA,       1,
  NA_character_,   27.9,   0,      5.3,    "mmol/L",  NA,       2,
  "Hyperglycemia", 27.9,   0,      5.3,    "mmol/L",  "4",      3,
  "Hyperglycemia", 27.8,   0,      5.3,    "mmol/L",  "3",      4,
  "Hyperglycemia", 14,     0,      5.3,    "mmol/L",  "3",      5,
  "Hyperglycemia", 13.9,   0,      5.3,    "mmol/L",  "0",      6,
  "Hyperglycemia", 5.3,    0,      NA,     "mmol/L",  "0",      7,
  # Unit missing cannot grade
  "Hyperglycemia", 5.3,    0,      5.3,    NA,        NA,       8,
  # AVAL missing cannot grade
  "Hyperglycemia", NA,     0,      5.3,    "mmol/L",  NA,       9,
)
input_glyci <- expected_glyci %>%
  select(-ATOXGRH)

expected_glyci_uscv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",     501,    0,       5.3,     "mg/dL",        NA,  1,
  NA_character_,    501,    0,       5.3,     "mg/dL",        NA,  2,
  "Hyperglycemia",  501,    0,       5.3,     "mg/dL",       "4",  3,
  "Hyperglycemia",  500,    0,       5.3,     "mg/dL",       "3",  4,
  "Hyperglycemia",  251,    0,       5.3,     "mg/dL",       "3",  5,
  "Hyperglycemia",  250,    0,       5.3,     "mg/dL",       "0",  6,
  "Hyperglycemia",  250,    0,       NA,      "mg/dL",       "0",  7,
  # Unit missing cannot grade
  "Hyperglycemia",  250,    0,       5.3,          NA,        NA,  8,
  # AVAL missing cannot grade
  "Hyperglycemia",  NA,     0,       5.3,     "mg/dL",        NA,  9,
)
input_glyci_uscv <- expected_glyci_uscv %>%
  select(-ATOXGRH)

## Test 56a: CTCAEv4 Hyperglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 56a: CTCAEv4 Hyperglycemia (SI unit)", {
  actual_glyci <- derive_var_atoxgr_dir(
    input_glyci,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glyci,
    compare = actual_glyci,
    keys = c("TESTNUM")
  )
})

## Test 56b: CTCAEv4 Hyperglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 56b: CTCAEv4 Hyperglycemia (USCV unit)", {
  actual_glyci_uscv <- derive_var_atoxgr_dir(
    input_glyci_uscv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glyci_uscv,
    compare = actual_glyci_uscv,
    keys = c("TESTNUM")
  )
})

## Test 56c: CTCAEv6 Hyperglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 56c: CTCAEv6 Hyperglycemia (SI unit)", {
  actual_glyci <- derive_var_atoxgr_dir(
    input_glyci,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glyci,
    compare = actual_glyci,
    keys = c("TESTNUM")
  )
})

## Test 56d: CTCAEv6 Hyperglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 56d: CTCAEv6 Hyperglycemia (USCV unit)", {
  actual_glyci_uscv <- derive_var_atoxgr_dir(
    input_glyci_uscv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glyci_uscv,
    compare = actual_glyci_uscv,
    keys = c("TESTNUM")
  )
})

### Hyperkalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI + CV unit
### Grade 4: >7.0 mmol/L
### Grade 3: >6.0 - 7.0 mmol/L
### Grade 2: >5.5 - 6.0 mmol/L
### Grade 1: >ULN - 5.5 mmol/L

expected_kalei <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    7.1,    0,      5.1,    "mmol/L",  NA,       1,
  NA_character_,   7.1,    0,      5.1,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hyperkalemia",  7.1,    0,      5.1,    "mmol/L",  "4",      3,
  "Hyperkalemia",  7,      0,      5.1,    "mmol/L",  "3",      4,
  "Hyperkalemia",  6.1,    0,      5.1,    "mmol/L",  "3",      5,
  "Hyperkalemia",  6,      0,      5.1,    "mmol/L",  "2",      6,
  "Hyperkalemia",  5.6,    0,      5.1,    "mmol/L",  "2",      7,
  "Hyperkalemia",  5.5,    0,      5.1,    "mmol/L",  "1",      8,
  "Hyperkalemia",  5.2,    0,      5.1,    "mmol/L",  "1",      9,
  "Hyperkalemia",  5.1,    0,      5.1,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hyperkalemia",  7.1,    0,      NA,     "mmol/L",  "4",      11,
  "Hyperkalemia",  7,      0,      NA,     "mmol/L",  "3",      12,
  "Hyperkalemia",  6.1,    0,      NA,     "mmol/L",  "3",      13,
  "Hyperkalemia",  6,      0,      NA,     "mmol/L",  "2",      14,
  "Hyperkalemia",  5.6,    0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hyperkalemia",  5.5,    0,      NA,     "mmol/L",  NA,       16,
  "Hyperkalemia",  5.2,    0,      NA,     "mmol/L",  NA,       17,
  "Hyperkalemia",  5.1,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hyperkalemia",  5.1,    0,      5.1,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hyperkalemia",  NA,     0,      5.1,    "mmol/L",  NA,       20,
)
input_kalei <- expected_kalei %>%
  select(-ATOXGRH)

## Test 57a: CTCAEv4 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 57a: CTCAEv4 Hyperkalemia (SI unit)", {
  actual_kalei <- derive_var_atoxgr_dir(
    input_kalei,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kalei,
    compare = actual_kalei,
    keys = c("TESTNUM")
  )
})

## Test 57b: CTCAEv4 Hyperkalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 57b: CTCAEv4 Hyperkalemia (USCV unit)", {
  actual_kalei <- derive_var_atoxgr_dir(
    input_kalei,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kalei,
    compare = actual_kalei,
    keys = c("TESTNUM")
  )
})

## Test 58a: CTCAEv5 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 58a: CTCAEv5 Hyperkalemia (SI unit)", {
  actual_kalei <- derive_var_atoxgr_dir(
    input_kalei,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kalei,
    compare = actual_kalei,
    keys = c("TESTNUM")
  )
})

## Test 58b: CTCAEv5 Hyperkalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 58b: CTCAEv5 Hyperkalemia (USCV unit)", {
  actual_kalei <- derive_var_atoxgr_dir(
    input_kalei,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kalei,
    compare = actual_kalei,
    keys = c("TESTNUM")
  )
})

## Test 58c: CTCAEv6 Hyperkalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 58c: CTCAEv6 Hyperkalemia (SI unit)", {
  actual_kalei <- derive_var_atoxgr_dir(
    input_kalei,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kalei,
    compare = actual_kalei,
    keys = c("TESTNUM")
  )
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

expected_magni_si <- tibble::tribble(
  ~ATOXDSCH,         ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",      3.4,   0,      0.8,    "mmol/L", NA,       1,
  NA_character_,     3.4,   0,      0.8,    "mmol/L", NA,       2,
  # ANRHI not missing
  "Hypermagnesemia", 3.4,   0,      0.8,    "mmol/L", "4",      3,
  "Hypermagnesemia", 3.3,   0,      0.8,    "mmol/L", "3",      4,
  "Hypermagnesemia", 1.24,  0,      0.8,    "mmol/L", "3",      5,
  "Hypermagnesemia", 1.23,  0,      0.8,    "mmol/L", "1",      6,
  "Hypermagnesemia", 0.81,  0,      0.8,    "mmol/L", "1",      7,
  "Hypermagnesemia", 0.8,   0,      0.8,    "mmol/L", "0",      8,
  # ANRHI missing - can grade 3-4
  "Hypermagnesemia", 3.4,   0,      NA,     "mmol/L", "4",      9,
  "Hypermagnesemia", 3.3,   0,      NA,     "mmol/L", "3",      10,
  "Hypermagnesemia", 1.24,  0,      NA,     "mmol/L", "3",      11,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypermagnesemia", 1.23,  0,      NA,     "mmol/L", NA,       12,
  "Hypermagnesemia", 0.81,  0,      NA,     "mmol/L", NA,       13,
  "Hypermagnesemia", 0.8,   0,      NA,     "mmol/L", NA,       14,
  # Unit missing cannot grade
  "Hypermagnesemia", 0.8,   0,      0.8,    NA,       NA,       15,
  # AVAL missing cannot grade
  "Hypermagnesemia", NA,    0,      0.8,    "mmol/L", NA,       16,
)
input_magni_si <- expected_magni_si %>%
  select(-ATOXGRH)

expected_magni_cv <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU, ~ATOXGRH,  ~TESTNUM,
  "Not a term",       8.1,    0,       1.2,     "mg/dL",       NA,  1,
  NA_character_,      8.1,    0,       1.2,     "mg/dL",       NA,  2,
  # ANRHI not missing
  "Hypermagnesemia",  8.1,    0,       1.2,     "mg/dL",      "4",  3,
  "Hypermagnesemia",  8,      0,       1.2,     "mg/dL",      "3",  4,
  "Hypermagnesemia",  3.1,    0,       1.2,     "mg/dL",      "3",  5,
  "Hypermagnesemia",  3,      0,       1.2,     "mg/dL",      "1",  6,
  "Hypermagnesemia",  1.21,   0,       1.2,     "mg/dL",      "1",  7,
  "Hypermagnesemia",  1.2,    0,       1.2,     "mg/dL",      "0",  8,
  # ANRHI missing - can grade 3-4
  "Hypermagnesemia",  8.1,    0,       NA,      "mg/dL",      "4",  9,
  "Hypermagnesemia",  8,      0,       NA,      "mg/dL",      "3",  10,
  "Hypermagnesemia",  3.1,    0,       NA,      "mg/dL",      "3",  11,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypermagnesemia",  3,      0,       NA,      "mg/dL",       NA,  12,
  "Hypermagnesemia",  1.21,   0,       NA,      "mg/dL",       NA,  13,
  "Hypermagnesemia",  1.2,    0,       NA,      "mg/dL",       NA,  14,
  # Unit missing cannot grade
  "Hypermagnesemia",  1.2,    0,       1.2,          NA,       NA,  15,
  # AVAL missing cannot grade
  "Hypermagnesemia",  NA,     0,       1.2,     "mg/dL",       NA,  16,
)
input_magni_cv <- expected_magni_cv %>%
  select(-ATOXGRH)


## Test 59a: CTCAEv4 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 59a: CTCAEv4 Hypermagnesemia (SI unit)", {
  actual_magni_si <- derive_var_atoxgr_dir(
    input_magni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_si,
    compare = actual_magni_si,
    keys = c("TESTNUM")
  )
})

## Test 59b: CTCAEv4 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 59b: CTCAEv4 Hypermagnesemia (USCV unit)", {
  actual_magni_cv <- derive_var_atoxgr_dir(
    input_magni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_cv,
    compare = actual_magni_cv,
    keys = c("TESTNUM")
  )
})

## Test 60a: CTCAEv5 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 60a: CTCAEv5 Hypermagnesemia (SI unit)", {
  actual_magni_si <- derive_var_atoxgr_dir(
    input_magni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_si,
    compare = actual_magni_si,
    keys = c("TESTNUM")
  )
})

## Test 60b: CTCAEv5 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 60b: CTCAEv5 Hypermagnesemia (USCV unit)", {
  actual_magni_cv <- derive_var_atoxgr_dir(
    input_magni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_cv,
    compare = actual_magni_cv,
    keys = c("TESTNUM")
  )
})

## Test 60c: CTCAEv6 Hypermagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 60c: CTCAEv6 Hypermagnesemia (SI unit)", {
  actual_magni_si <- derive_var_atoxgr_dir(
    input_magni_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_si,
    compare = actual_magni_si,
    keys = c("TESTNUM")
  )
})

## Test 60d: CTCAEv6 Hypermagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 60d: CTCAEv6 Hypermagnesemia (USCV unit)", {
  actual_magni_cv <- derive_var_atoxgr_dir(
    input_magni_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni_cv,
    compare = actual_magni_cv,
    keys = c("TESTNUM")
  )
})

### Hypernatremia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI + CV unit
### Grade 4: >160 mmol/L
### Grade 3: >155 - 160 mmol/L
### Grade 2: >150 - 155 mmol/L
### Grade 1: >ULN - 150 mmol/L

expected_natri <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",     161,    0,      140,    "mmol/L",  NA,       1,
  NA_character_,    161,    0,      140,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypernatremia",  161,    0,      140,    "mmol/L",  "4",      3,
  "Hypernatremia",  160,    0,      140,    "mmol/L",  "3",      4,
  "Hypernatremia",  156,    0,      140,    "mmol/L",  "3",      5,
  "Hypernatremia",  155,    0,      140,    "mmol/L",  "2",      6,
  "Hypernatremia",  151,    0,      140,    "mmol/L",  "2",      7,
  "Hypernatremia",  150,    0,      140,    "mmol/L",  "1",      8,
  "Hypernatremia",  141,    0,      140,    "mmol/L",  "1",      9,
  "Hypernatremia",  140,    0,      140,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 3-4
  "Hypernatremia",  161,    0,      NA,     "mmol/L",  "4",      11,
  "Hypernatremia",  160,    0,      NA,     "mmol/L",  "3",      12,
  "Hypernatremia",  156,    0,      NA,     "mmol/L",  "3",      13,
  "Hypernatremia",  155,    0,      NA,     "mmol/L",  "2",      14,
  "Hypernatremia",  151,    0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypernatremia",  150,    0,      NA,     "mmol/L",  NA,       16,
  "Hypernatremia",  141,    0,      NA,     "mmol/L",  NA,       17,
  "Hypernatremia",  140,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypernatremia",  140,    0,      140,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypernatremia",  NA,     0,      140,    "mmol/L",  NA,       20,
)
input_natri <- expected_natri %>%
  select(-ATOXGRH)

## Test 61a: CTCAEv4 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 61a: CTCAEv4 Hypernatremia (SI unit)", {
  actual_natri <- derive_var_atoxgr_dir(
    input_natri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natri,
    compare = actual_natri,
    keys = c("TESTNUM")
  )
})

## Test 61b: CTCAEv4 Hypernatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 61b: CTCAEv4 Hypernatremia (USCV unit)", {
  actual_natri <- derive_var_atoxgr_dir(
    input_natri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natri,
    compare = actual_natri,
    keys = c("TESTNUM")
  )
})

## Test 62a: CTCAEv5 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 62a: CTCAEv5 Hypernatremia (SI unit)", {
  actual_natri <- derive_var_atoxgr_dir(
    input_natri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natri,
    compare = actual_natri,
    keys = c("TESTNUM")
  )
})

## Test 62b: CTCAEv5 Hypernatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 62b: CTCAEv5 Hypernatremia (USCV unit)", {
  actual_natri <- derive_var_atoxgr_dir(
    input_natri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natri,
    compare = actual_natri,
    keys = c("TESTNUM")
  )
})

## Test 62c: CTCAEv6 Hypernatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 62c: CTCAEv6 Hypernatremia (SI unit)", {
  actual_natri <- derive_var_atoxgr_dir(
    input_natri,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natri,
    compare = actual_natri,
    keys = c("TESTNUM")
  )
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

expected_trigi_si <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
  "Not a term",            11.5,   0,      2.1,    "mmol/L",  NA,
  NA_character_,           11.5,   0,      2.1,    "mmol/L",  NA,
  "Hypertriglyceridemia",  11.5,   0,      2.1,    "mmol/L",  "4",
  "Hypertriglyceridemia",  11.4,   0,      2.1,    "mmol/L",  "3",
  "Hypertriglyceridemia",  5.8,    0,      2.1,    "mmol/L",  "3",
  "Hypertriglyceridemia",  5.7,    0,      2.1,    "mmol/L",  "2",
  "Hypertriglyceridemia",  3.43,   0,      2.1,    "mmol/L",  "2",
  "Hypertriglyceridemia",  3.42,   0,      2.1,    "mmol/L",  "1",
  "Hypertriglyceridemia",  1.72,   0,      2.1,    "mmol/L",  "1",
  "Hypertriglyceridemia",  1.71,   0,      2.1,    "mmol/L",  "0",
  # Unit missing cannot grade
  "Hypertriglyceridemia",  1.71,   0,      2.1,    NA,        NA,
  # AVAL missing cannot grade
  "Hypertriglyceridemia",  NA,     0,      2.1,    "mmol/L",  NA,
)
input_trigi_si <- expected_trigi_si %>%
  select(-ATOXGRH)

expected_trigi_cv <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~ANRLO, ~ANRHI,  ~AVALU,   ~ATOXGRH,
  "Not a term",            1001,   0,      50,      "mg/dL",  NA,
  NA_character_,           1001,   0,      50,      "mg/dL",  NA,
  "Hypertriglyceridemia",  1001,   0,      50,      "mg/dL",  "4",
  "Hypertriglyceridemia",  1000,   0,      50,      "mg/dL",  "3",
  "Hypertriglyceridemia",  501,    0,      50,      "mg/dL",  "3",
  "Hypertriglyceridemia",  500,    0,      50,      "mg/dL",  "2",
  "Hypertriglyceridemia",  301,    0,      50,      "mg/dL",  "2",
  "Hypertriglyceridemia",  300,    0,      50,      "mg/dL",  "1",
  "Hypertriglyceridemia",  151,    0,      50,      "mg/dL",  "1",
  "Hypertriglyceridemia",  150,    0,      50,      "mg/dL",  "0",
  # Unit missing cannot grade
  "Hypertriglyceridemia",  150,    0,      50,      NA,       NA,
  # AVAL missing cannot grade
  "Hypertriglyceridemia",  NA,     0,      50,      "mg/dL",  NA,
)
input_trigi_cv <- expected_trigi_cv %>%
  select(-ATOXGRH)


## Test 63a: CTCAEv4 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 63a: CTCAEv4 Hypertriglyceridemia (SI unit)", {
  actual_trigi_si <- derive_var_atoxgr_dir(
    input_trigi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_si,
    compare = actual_trigi_si,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 63b: CTCAEv4 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 63b: CTCAEv4 Hypertriglyceridemia (USCV unit)", {
  actual_trigi_cv <- derive_var_atoxgr_dir(
    input_trigi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_cv,
    compare = actual_trigi_cv,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 64a: CTCAEv5 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 64a: CTCAEv5 Hypertriglyceridemia (SI unit)", {
  actual_trigi_si <- derive_var_atoxgr_dir(
    input_trigi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_si,
    compare = actual_trigi_si,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 64b: CTCAEv5 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 64b: CTCAEv5 Hypertriglyceridemia (USCV unit)", {
  actual_trigi_cv <- derive_var_atoxgr_dir(
    input_trigi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_cv,
    compare = actual_trigi_cv,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 64c: CTCAEv6 Hypertriglyceridemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 64c: CTCAEv6 Hypertriglyceridemia (SI unit)", {
  actual_trigi_si <- derive_var_atoxgr_dir(
    input_trigi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_si,
    compare = actual_trigi_si,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 64d: CTCAEv6 Hypertriglyceridemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 64d: CTCAEv6 Hypertriglyceridemia (USCV unit)", {
  actual_trigi_cv <- derive_var_atoxgr_dir(
    input_trigi_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi_cv,
    compare = actual_trigi_cv,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

### Hyperuricemia (NCICTCAEv4)
### NCICTCAEv5 only has grade 3
### SI unit is umol/L = 1000 * mmol/L
### Grade 4: >0.59 mmol/L;
### Grade 3: >ULN - 0.59 mmol/L

### CV unit is mg/dL
### Grade 4: >10 mg/dL
### Grade 3: >ULN - 10 mg/dL


expected_urici_si <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,    ~AVALU, ~ATOXGRH,
  "Not a term",     591,    0,       200,     "umol/L",       NA,
  NA_character_,    591,    0,       200,     "umol/L",       NA,
  # ANRHI not missing
  "Hyperuricemia",  591,    0,       200,     "umol/L",      "4",
  "Hyperuricemia",  590,    0,       200,     "umol/L",      "3",
  "Hyperuricemia",  201,    0,       200,     "umol/L",      "3",
  "Hyperuricemia",  200,    0,       200,     "umol/L",      "0",
  # ANRHI missing - can grade 4
  "Hyperuricemia",  591,    0,       NA,      "umol/L",      "4",
  # ANRHI missing - can NOT grade 0 or 3
  "Hyperuricemia",  590,    0,       NA,      "umol/L",       NA,
  "Hyperuricemia",  201,    0,       NA,      "umol/L",       NA,
  "Hyperuricemia",  200,    0,       NA,      "umol/L",       NA,
  # Unit missing cannot grade
  "Hyperuricemia",  200,    0,       200,           NA,      "0",
  # AVAL missing cannot grade
  "Hyperuricemia",  NA,     0,       200,     "umol/L",       NA,
)
input_urici_si <- expected_urici_si %>%
  select(-ATOXGRH)

expected_urici_cv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL, ~ANRLO, ~ANRHI,  ~AVALU,  ~ATOXGRH,
  "Not a term",     11,    0,      5,      "mg/dL",  NA,
  NA_character_,    11,    0,      5,      "mg/dL",  NA,
  # ANRHI not missing
  "Hyperuricemia",  11,    0,      5,      "mg/dL",  "4",
  "Hyperuricemia",  10,    0,      5,      "mg/dL",  "3",
  "Hyperuricemia",  6,     0,      5,      "mg/dL",  "3",
  "Hyperuricemia",  5,     0,      5,      "mg/dL",  "0",
  # ANRHI missing - can grade 4
  "Hyperuricemia",  11,    0,      NA,     "mg/dL",  "4",
  # ANRHI missing - can NOT grade 0 or 3
  "Hyperuricemia",  10,    0,      NA,     "mg/dL",  NA,
  "Hyperuricemia",  6,     0,      NA,     "mg/dL",  NA,
  "Hyperuricemia",  5,     0,      NA,     "mg/dL",  NA,
  # Unit missing cannot grade
  "Hyperuricemia",  5,     0,      5,           NA,  "0",
  # AVAL missing cannot grade
  "Hyperuricemia",  NA,    0,      5,      "mg/dL",  NA,
)

input_urici_cv <- expected_urici_cv %>%
  select(-ATOXGRH)

### Hyperuricemia (NCICTCAEv5)
### NCICTCAEv5 only has grade 3
### Grade 3: >ULN

## Test 65a: CTCAEv5 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 65a: CTCAEv5 Hyperuricemia (SI unit)", {
  expected_urici_si <- expected_urici_si %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")
  input_urici_si <- expected_urici_si %>%
    select(-ATOXGRH)

  actual_urici_si <- derive_var_atoxgr_dir(
    input_urici_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_si,
    compare = actual_urici_si,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 65b: CTCAEv5 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 65b: CTCAEv5 Hyperuricemia (USCV unit)", {
  expected_urici_cv <- expected_urici_cv %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")
  input_urici_cv <- expected_urici_cv %>%
    select(-ATOXGRH)

  actual_urici_cv <- derive_var_atoxgr_dir(
    input_urici_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_cv,
    compare = actual_urici_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 65c: CTCAEv6 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 65c: CTCAEv6 Hyperuricemia (SI unit)", {
  expected_urici_si <- expected_urici_si %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")
  input_urici_si <- expected_urici_si %>%
    select(-ATOXGRH)

  actual_urici_si <- derive_var_atoxgr_dir(
    input_urici_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_si,
    compare = actual_urici_si,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 65d: CTCAEv6 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 65d: CTCAEv6 Hyperuricemia (USCV unit)", {
  expected_urici_cv <- expected_urici_cv %>%
    filter(is.na(ATOXGRH) | ATOXGRH != "4")
  input_urici_cv <- expected_urici_cv %>%
    select(-ATOXGRH)

  actual_urici_cv <- derive_var_atoxgr_dir(
    input_urici_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_cv,
    compare = actual_urici_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})


# If unit missing then grade CANNOT be calculated as needed for grade 4
## Test 66a: CTCAEv4 Hyperuricemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 66a: CTCAEv4 Hyperuricemia (SI unit)", {
  expected_urici_si <- expected_urici_si %>%
    mutate(ATOXGRH = if_else(is.na(AVALU), NA_character_, ATOXGRH))

  input_urici_si <- expected_urici_si %>%
    select(-ATOXGRH)

  actual_urici_si <- derive_var_atoxgr_dir(
    input_urici_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_si,
    compare = actual_urici_si,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

## Test 66b: CTCAEv4 Hyperuricemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 66b: CTCAEv4 Hyperuricemia (USCV unit)", {
  expected_urici_cv <- expected_urici_cv %>%
    mutate(ATOXGRH = if_else(is.na(AVALU), NA_character_, ATOXGRH))

  input_urici_cv <- expected_urici_cv %>%
    select(-ATOXGRH)

  actual_urici_cv <- derive_var_atoxgr_dir(
    input_urici_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_cv,
    compare = actual_urici_cv,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### Hypoalbuminemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI init is g/L - CV unit is g/dL
### Grade 3: <20 g/L
### Grade 2: <30 - 20 g/L
### Grade 1: <LLN - 30 g/L

expected_albd_si <- tibble::tribble(
  ~ATOXDSCL,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL,
  "Not a term",       19,     40,     100,    "G/L",  NA,
  NA_character_,      19,     40,     100,    "G/L",  NA,
  # ANRLO not missing
  "Hypoalbuminemia",  19,     40,     100,    "G/L",  "3",
  "Hypoalbuminemia",  20,     40,     100,    "G/L",  "2",
  "Hypoalbuminemia",  29,     40,     100,    "G/L",  "2",
  "Hypoalbuminemia",  30,     40,     100,    "G/L",  "1",
  "Hypoalbuminemia",  39,     40,     100,    "G/L",  "1",
  "Hypoalbuminemia",  40,     40,     100,    "G/L",  "0",
  "Hypoalbuminemia",  40,     40,     NA,     "G/L",  "0",
  # ANRLO missing - can grade 2-3
  "Hypoalbuminemia",  19,     NA,     100,    "G/L",  "3",
  "Hypoalbuminemia",  20,     NA,     100,    "G/L",  "2",
  "Hypoalbuminemia",  29,     NA,     100,    "G/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypoalbuminemia",  30,     NA,     100,    "G/L",  NA,
  "Hypoalbuminemia",  39,     NA,     100,    "G/L",  NA,
  "Hypoalbuminemia",  40,     NA,     100,    "G/L",  NA,
  # Unit missing cannot grade
  "Hypoalbuminemia",  40,     40,     100,    NA,     NA,
  # AVAL missing cannot grade
  "Hypoalbuminemia",  NA,     40,     100,    "G/L",  NA,
)
input_albd_si <- expected_albd_si %>%
  select(-ATOXGRL)

expected_albd_cv <- expected_albd_si %>%
  mutate(
    AVAL = AVAL / 10,
    ANRLO = ANRLO / 10,
    ANRHI = ANRHI / 10,
    AVALU = if_else(is.na(AVALU), NA_character_, "g/dL")
  )

input_albd_cv <- expected_albd_cv %>%
  select(-ATOXGRL)


## Test 67a: CTCAEv4 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 67a: CTCAEv4 Hypoalbuminemia (SI unit)", {
  actual_albd_si <- derive_var_atoxgr_dir(
    input_albd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_si,
    compare = actual_albd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 67b: CTCAEv4 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 67b: CTCAEv4 Hypoalbuminemia (USCV unit)", {
  actual_albd_cv <- derive_var_atoxgr_dir(
    input_albd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_cv,
    compare = actual_albd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 68a: CTCAEv5 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 68a: CTCAEv5 Hypoalbuminemia (SI unit)", {
  actual_albd_si <- derive_var_atoxgr_dir(
    input_albd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_si,
    compare = actual_albd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 68b: CTCAEv5 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 68b: CTCAEv5 Hypoalbuminemia (USCV unit)", {
  actual_albd_cv <- derive_var_atoxgr_dir(
    input_albd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_cv,
    compare = actual_albd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 68c: CTCAEv6 Hypoalbuminemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 68c: CTCAEv6 Hypoalbuminemia (SI unit)", {
  actual_albd_si <- derive_var_atoxgr_dir(
    input_albd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_si,
    compare = actual_albd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 68d: CTCAEv6 Hypoalbuminemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 68d: CTCAEv6 Hypoalbuminemia (USCV unit)", {
  actual_albd_cv <- derive_var_atoxgr_dir(
    input_albd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd_cv,
    compare = actual_albd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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

expected_calcd_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",    1.4,    4,      100,    "mmol/L",  NA,
  NA_character_,   1.4,    4,      100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hypocalcemia",  1.4,    4,      100,    "mmol/L",  "4",
  "Hypocalcemia",  1.5,    4,      100,    "mmol/L",  "3",
  "Hypocalcemia",  1.74,   4,      100,    "mmol/L",  "3",
  "Hypocalcemia",  1.75,   4,      100,    "mmol/L",  "2",
  "Hypocalcemia",  1.9,    4,      100,    "mmol/L",  "2",
  "Hypocalcemia",  2,      4,      100,    "mmol/L",  "1",
  "Hypocalcemia",  3,      4,      100,    "mmol/L",  "1",
  "Hypocalcemia",  4,      4,      100,    "mmol/L",  "0",
  # ANRLO missing - can grade 2-4
  "Hypocalcemia",  1.4,    4,      NA,     "mmol/L",  "4",
  "Hypocalcemia",  1.5,    4,      NA,     "mmol/L",  "3",
  "Hypocalcemia",  1.74,   4,      NA,     "mmol/L",  "3",
  "Hypocalcemia",  1.75,   4,      NA,     "mmol/L",  "2",
  "Hypocalcemia",  1.9,    4,      NA,     "mmol/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia",  2,      4,      NA,     "mmol/L",  "1",
  "Hypocalcemia",  3,      4,      NA,     "mmol/L",  "1",
  "Hypocalcemia",  4,      4,      NA,     "mmol/L",  "0",
  # Unit missing cannot grade
  "Hypocalcemia",  4,      4,      100,    NA,        NA,
  # AVAL missing cannot grade
  "Hypocalcemia",  NA,     4,      100,    "mmol/L",  NA,
)
input_calcd_si <- expected_calcd_si %>%
  select(-ATOXGRL)

expected_calcd_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL,
  "Not a term",    5.9,    9,       100,     "mg/dL",        NA,
  NA_character_,   5.9,    9,       100,     "mg/dL",        NA,
  # ANRLO not missing
  "Hypocalcemia",  5.9,    9,       100,     "mg/dL",       "4",
  "Hypocalcemia",  6,      9,       100,     "mg/dL",       "3",
  "Hypocalcemia",  6.9,    9,       100,     "mg/dL",       "3",
  "Hypocalcemia",  7,      9,       100,     "mg/dL",       "2",
  "Hypocalcemia",  7.9,    9,       100,     "mg/dL",       "2",
  "Hypocalcemia",  8,      9,       100,     "mg/dL",       "1",
  "Hypocalcemia",  8.9,    9,       100,     "mg/dL",       "1",
  "Hypocalcemia",  9,      9,       100,     "mg/dL",       "0",
  # ANRLO missing - can grade 2-4
  "Hypocalcemia",  5.9,    9,       NA,      "mg/dL",       "4",
  "Hypocalcemia",  6,      9,       NA,      "mg/dL",       "3",
  "Hypocalcemia",  6.9,    9,       NA,      "mg/dL",       "3",
  "Hypocalcemia",  7,      9,       NA,      "mg/dL",       "2",
  "Hypocalcemia",  7.9,    9,       NA,      "mg/dL",       "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia",  8,      9,       NA,      "mg/dL",       "1",
  "Hypocalcemia",  8.9,    9,       NA,      "mg/dL",       "1",
  "Hypocalcemia",  9,      9,       NA,      "mg/dL",       "0",
  # Unit missing cannot grade
  "Hypocalcemia",  9,      9,       100,          NA,        NA,
  # AVAL missing cannot grade
  "Hypocalcemia",  NA,     9,       100,     "mg/dL",        NA,
)
input_calcd_cv <- expected_calcd_cv %>%
  select(-ATOXGRL)

## Test 69a: CTCAEv4 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 69a: CTCAEv4 Hypocalcemia (SI unit)", {
  actual_calcd_si <- derive_var_atoxgr_dir(
    input_calcd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_si,
    compare = actual_calcd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 69b: CTCAEv4 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 69b: CTCAEv4 Hypocalcemia (USCV unit)", {
  actual_calcd_cv <- derive_var_atoxgr_dir(
    input_calcd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_cv,
    compare = actual_calcd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 70a: CTCAEv5 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 70a: CTCAEv5 Hypocalcemia (SI unit)", {
  actual_calcd_si <- derive_var_atoxgr_dir(
    input_calcd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_si,
    compare = actual_calcd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 70b: CTCAEv5 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 70b: CTCAEv5 Hypocalcemia (USCV unit)", {
  actual_calcd_cv <- derive_var_atoxgr_dir(
    input_calcd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_cv,
    compare = actual_calcd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 70c: CTCAEv6 Hypocalcemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 70c: CTCAEv6 Hypocalcemia (SI unit)", {
  actual_calcd_si <- derive_var_atoxgr_dir(
    input_calcd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_si,
    compare = actual_calcd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 70d: CTCAEv6 Hypocalcemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 70d: CTCAEv6 Hypocalcemia (USCV unit)", {
  actual_calcd_cv <- derive_var_atoxgr_dir(
    input_calcd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_cv,
    compare = actual_calcd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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


expected_caliond_si <- tibble::tribble(
  ~ATOXDSCL,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",              0.79,   1.3,    100,    "mmol/L",  NA,
  NA_character_,             0.79,   1.3,    100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hypocalcemia (Ionized)",  0.79,   1.3,    100,    "mmol/L",  "4",
  "Hypocalcemia (Ionized)",  0.8,    1.3,    100,    "mmol/L",  "3",
  "Hypocalcemia (Ionized)",  0.89,   1.3,    100,    "mmol/L",  "3",
  "Hypocalcemia (Ionized)",  0.9,    1.3,    100,    "mmol/L",  "2",
  "Hypocalcemia (Ionized)",  0.99,   1.3,    100,    "mmol/L",  "2",
  "Hypocalcemia (Ionized)",  1,      1.3,    100,    "mmol/L",  "1",
  "Hypocalcemia (Ionized)",  1.29,   1.3,    100,    "mmol/L",  "1",
  "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    "mmol/L",  "0",
  # ANRLO missing - can grade 2-4
  "Hypocalcemia (Ionized)",  0.79,   NA,     100,    "mmol/L",  "4",
  "Hypocalcemia (Ionized)",  0.8,    NA,     100,    "mmol/L",  "3",
  "Hypocalcemia (Ionized)",  0.89,   NA,     100,    "mmol/L",  "3",
  "Hypocalcemia (Ionized)",  0.9,    NA,     100,    "mmol/L",  "2",
  "Hypocalcemia (Ionized)",  0.99,   NA,     100,    "mmol/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia (Ionized)",  1,      1.3,    NA,     "mmol/L",  "1",
  "Hypocalcemia (Ionized)",  1.29,   1.3,    NA,     "mmol/L",  "1",
  "Hypocalcemia (Ionized)",  1.3,    1.3,    NA,     "mmol/L",  "0",
  # Unit missing cannot grade
  "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    NA,        NA,
  # AVAL missing cannot grade
  "Hypocalcemia (Ionized)",  NA,     1.3,    100,    "mmol/L",  NA,
)
input_caliond_si <- expected_caliond_si %>%
  select(-ATOXGRL)

expected_caliond_cv <- expected_caliond_si %>%
  mutate(
    AVAL = AVAL * 4,
    ANRHI = ANRHI * 4,
    ANRLO = ANRLO * 4,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )
input_caliond_cv <- expected_caliond_cv %>%
  select(-ATOXGRL)

## Test 71a: CTCAEv4 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 71a: CTCAEv4 Hypocalcemia (Ionized) (SI unit)", {
  actual_caliond_si <- derive_var_atoxgr_dir(
    input_caliond_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_si,
    compare = actual_caliond_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 71b: CTCAEv4 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 71b: CTCAEv4 Hypocalcemia (Ionized) (USCV unit)", {
  actual_caliond_cv <- derive_var_atoxgr_dir(
    input_caliond_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_cv,
    compare = actual_caliond_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 72a: CTCAEv5 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 72a: CTCAEv5 Hypocalcemia (Ionized) (SI unit)", {
  actual_caliond_si <- derive_var_atoxgr_dir(
    input_caliond_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_si,
    compare = actual_caliond_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 72b: CTCAEv5 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 72b: CTCAEv5 Hypocalcemia (Ionized) (USCV unit)", {
  actual_caliond_cv <- derive_var_atoxgr_dir(
    input_caliond_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_cv,
    compare = actual_caliond_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 72c: CTCAEv6 Hypocalcemia (Ionized) (SI unit) ----
test_that("derive_var_atoxgr_dir Test 72c: CTCAEv6 Hypocalcemia (Ionized) (SI unit)", {
  actual_caliond_si <- derive_var_atoxgr_dir(
    input_caliond_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_si,
    compare = actual_caliond_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 72d: CTCAEv6 Hypocalcemia (Ionized) (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 72d: CTCAEv6 Hypocalcemia (Ionized) (USCV unit)", {
  actual_caliond_cv <- derive_var_atoxgr_dir(
    input_caliond_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_cv,
    compare = actual_caliond_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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

expected_glycd_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,  ~ANRHI,    ~AVALU,  ~ATOXGRL,
  "Not a term",    1.69,   4,       100,     "mmol/L",        NA,
  NA_character_,   1.69,   4,       100,     "mmol/L",        NA,
  # ANRLO not missing
  "Hypoglycemia",  1.69,   4,       100,     "mmol/L",       "4",
  "Hypoglycemia",  1.7,    4,       100,     "mmol/L",       "3",
  "Hypoglycemia",  2.19,   4,       100,     "mmol/L",       "3",
  "Hypoglycemia",  2.2,    4,       100,     "mmol/L",       "2",
  "Hypoglycemia",  2.9,    4,       100,     "mmol/L",       "2",
  "Hypoglycemia",  3,      4,       100,     "mmol/L",       "1",
  "Hypoglycemia",  3.9,    4,       100,     "mmol/L",       "1",
  "Hypoglycemia",  4,      4,       100,     "mmol/L",       "0",
  # ANRLO missing - can grade 2-4
  "Hypoglycemia",  1.69,   NA,      100,     "mmol/L",       "4",
  "Hypoglycemia",  1.7,    NA,      100,     "mmol/L",       "3",
  "Hypoglycemia",  2.19,   NA,      100,     "mmol/L",       "3",
  "Hypoglycemia",  2.2,    NA,      100,     "mmol/L",       "2",
  "Hypoglycemia",  2.9,    NA,      100,     "mmol/L",       "2",
  # ANRLO missing - can NOT grade 0 or 1mg/dL
  "Hypoglycemia",  3,      NA,      100,     "mmol/L",        NA,
  "Hypoglycemia",  3.9,    NA,      100,     "mmol/L",        NA,
  "Hypoglycemia",  4,      NA,      100,     "mmol/L",        NA,
  # Unit missing cannot grade
  "Hypoglycemia",  4,      4,       100,           NA,        NA,
  # AVAL missing cannot grade
  "Hypoglycemia",  NA,     4,       100,     "mmol/L",        NA,
)
input_glycd_si <- expected_glycd_si %>%
  select(-ATOXGRL)

expected_glycd_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI,   ~AVALU,  ~ATOXGRL,
  "Not a term",    29,     70,     100,     "mg/dL",        NA,
  NA_character_,   29,     70,     100,      "g/dL",        NA,
  # ANRLO not missing
  "Hypoglycemia",  29,     70,     100,     "mg/dL",       "4",
  "Hypoglycemia",  30,     70,     100,     "mg/dL",       "3",
  "Hypoglycemia",  39,     70,     100,     "mg/dL",       "3",
  "Hypoglycemia",  40,     70,     100,     "mg/dL",       "2",
  "Hypoglycemia",  54,     70,     100,     "mg/dL",       "2",
  "Hypoglycemia",  55,     70,     100,     "mg/dL",       "1",
  "Hypoglycemia",  69,     70,     100,     "mg/dL",       "1",
  "Hypoglycemia",  70,     70,     100,     "mg/dL",       "0",
  # ANRLO missing - can grade 2-4
  "Hypoglycemia",  29,     NA,     100,     "mg/dL",       "4",
  "Hypoglycemia",  30,     NA,     100,     "mg/dL",       "3",
  "Hypoglycemia",  39,     NA,     100,     "mg/dL",       "3",
  "Hypoglycemia",  40,     NA,     100,     "mg/dL",       "2",
  "Hypoglycemia",  54,     NA,     100,     "mg/dL",       "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypoglycemia",  55,     NA,     100,     "mg/dL",        NA,
  "Hypoglycemia",  69,     NA,     100,     "mg/dL",        NA,
  "Hypoglycemia",  70,     NA,     100,     "mg/dL",        NA,
  # Unit missing cannot     grade
  "Hypoglycemia",  70,     70,     100,          NA,        NA,
  # AVAL missing cannot grade
  "Hypoglycemia",  NA,     70,     100,     "mg/dL",        NA,
)
input_glycd_cv <- expected_glycd_cv %>%
  select(-ATOXGRL)

## Test 73a: CTCAEv4 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 73a: CTCAEv4 Hypoglycemia (SI unit)", {
  actual_glycd_si <- derive_var_atoxgr_dir(
    input_glycd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_si,
    compare = actual_glycd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 73b: CTCAEv4 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 73b: CTCAEv4 Hypoglycemia (USCV unit)", {
  actual_glycd_cv <- derive_var_atoxgr_dir(
    input_glycd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_cv,
    compare = actual_glycd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})


## Test 74a: CTCAEv5 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 74a: CTCAEv5 Hypoglycemia (SI unit)", {
  actual_glycd_si <- derive_var_atoxgr_dir(
    input_glycd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_si,
    compare = actual_glycd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 74b: CTCAEv5 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 74b: CTCAEv5 Hypoglycemia (USCV unit)", {
  actual_glycd_cv <- derive_var_atoxgr_dir(
    input_glycd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_cv,
    compare = actual_glycd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 74c: CTCAEv6 Hypoglycemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 74c: CTCAEv6 Hypoglycemia (SI unit)", {
  actual_glycd_si <- derive_var_atoxgr_dir(
    input_glycd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_si,
    compare = actual_glycd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 74d: CTCAEv6 Hypoglycemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 74d: CTCAEv6 Hypoglycemia (USCV unit)", {
  actual_glycd_cv <- derive_var_atoxgr_dir(
    input_glycd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd_cv,
    compare = actual_glycd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hypokalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### SI unit and CV unit is mmol/L
### Grade 4: <2.5 mmol/L
### Grade 3: <3.0 - 2.5 mmol/L
### Grade 2: <LLN - 3.0 mmol/L

expected_kaled <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",   2.49,   4,      100,    "mmol/L",  NA,
  NA_character_,  2.49,   4,      100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hypokalemia",  2.49,   4,      100,    "mmol/L",  "4",
  "Hypokalemia",  2.5,    4,      100,    "mmol/L",  "3",
  "Hypokalemia",  2.9,    4,      100,    "mmol/L",  "3",
  "Hypokalemia",  3,      4,      100,    "mmol/L",  "2",
  "Hypokalemia",  3.9,    4,      100,    "mmol/L",  "2",
  "Hypokalemia",  4,      4,      100,    "mmol/L",  "0",
  # ANRLO missing - can grade 3-4
  "Hypokalemia",  2.49,   NA,     100,    "mmol/L",  "4",
  "Hypokalemia",  2.5,    NA,     100,    "mmol/L",  "3",
  "Hypokalemia",  2.9,    NA,     100,    "mmol/L",  "3",
  # ANRLO missing - can NOT grade 0 or 2
  "Hypokalemia",  3,      NA,     100,    "mmol/L",  NA,
  "Hypokalemia",  3.9,    NA,     100,    "mmol/L",  NA,
  "Hypokalemia",  4,      NA,     NA,     "mmol/L",  NA,
  # Unit missing cannot grade
  "Hypokalemia",  4,      4,      100,    NA,        NA,
  # AVAL missing cannot grade
  "Hypokalemia",  NA,     4,      100,    "mmol/L",  NA,
)
input_kaled <- expected_kaled %>%
  select(-ATOXGRL)

## Test 75a: CTCAEv4 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 75a: CTCAEv4 Hypokalemia (SI unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 75b: CTCAEv4 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 75b: CTCAEv4 Hypokalemia (USCV unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 76a: CTCAEv5 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 76a: CTCAEv5 Hypokalemia (SI unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 76b: CTCAEv5 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 76b: CTCAEv5 Hypokalemia (USCV unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 76c: CTCAEv6 Hypokalemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 76c: CTCAEv6 Hypokalemia (SI unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 76d: CTCAEv6 Hypokalemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 76d: CTCAEv6 Hypokalemia (USCV unit)", {
  actual_kaled <- derive_var_atoxgr_dir(
    input_kaled,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_kaled,
    compare = actual_kaled,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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

expected_magnd_si <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",      0.29,   1,      100,    "mmol/L",  NA,
  NA_character_,     0.29,   1,      100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hypomagnesemia",  0.29,   1,      100,    "mmol/L",  "4",
  "Hypomagnesemia",  0.3,    1,      100,    "mmol/L",  "3",
  "Hypomagnesemia",  0.39,   1,      100,    "mmol/L",  "3",
  "Hypomagnesemia",  0.4,    1,      100,    "mmol/L",  "2",
  "Hypomagnesemia",  0.49,   1,      100,    "mmol/L",  "2",
  "Hypomagnesemia",  0.5,    1,      100,    "mmol/L",  "1",
  "Hypomagnesemia",  0.9,    1,      100,    "mmol/L",  "1",
  "Hypomagnesemia",  1,      1,      100,    "mmol/L",  "0",
  # ANRLO missing - can grade 2-4
  "Hypomagnesemia",  0.29,   NA,     100,    "mmol/L",  "4",
  "Hypomagnesemia",  0.3,    NA,     100,    "mmol/L",  "3",
  "Hypomagnesemia",  0.39,   NA,     100,    "mmol/L",  "3",
  "Hypomagnesemia",  0.4,    NA,     100,    "mmol/L",  "2",
  "Hypomagnesemia",  0.49,   NA,     100,    "mmol/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypomagnesemia",  0.5,    NA,     100,    "mmol/L",  NA,
  "Hypomagnesemia",  0.9,    NA,     100,    "mmol/L",  NA,
  "Hypomagnesemia",  1,      NA,     100,    "mmol/L",  NA,
  # Unit missing cannot grade
  "Hypomagnesemia",  1,      1,      100,    NA,        NA,
  # AVAL missing cannot grade
  "Hypomagnesemia",  NA,     1,      100,    "mmol/L",  NA,
)
input_magnd_si <- expected_magnd_si %>%
  select(-ATOXGRL)

expected_magnd_cv <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL,
  "Not a term",      0.69,   1.5,     10,      "mg/dL",        NA,
  NA_character_,     0.69,   1.5,     10,      "mg/dL",        NA,
  # ANRLO not missing
  "Hypomagnesemia",  0.69,   1.5,     10,      "mg/dL",       "4",
  "Hypomagnesemia",  0.7,    1.5,     10,      "mg/dL",       "3",
  "Hypomagnesemia",  0.89,   1.5,     10,      "mg/dL",       "3",
  "Hypomagnesemia",  0.9,    1.5,     10,      "mg/dL",       "2",
  "Hypomagnesemia",  1.19,   1.5,     10,      "mg/dL",       "2",
  "Hypomagnesemia",  1.2,    1.5,     10,      "mg/dL",       "1",
  "Hypomagnesemia",  1.49,   1.5,     10,      "mg/dL",       "1",
  "Hypomagnesemia",  1.5,    1.5,     10,      "mg/dL",       "0",
  # ANRLO missing - can grade 2-4
  "Hypomagnesemia",  0.69,   NA,      10,      "mg/dL",       "4",
  "Hypomagnesemia",  0.7,    NA,      10,      "mg/dL",       "3",
  "Hypomagnesemia",  0.89,   NA,      10,      "mg/dL",       "3",
  "Hypomagnesemia",  0.9,    NA,      10,      "mg/dL",       "2",
  "Hypomagnesemia",  1.19,   NA,      10,      "mg/dL",       "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypomagnesemia",  1.2,    NA,      10,      "mg/dL",        NA,
  "Hypomagnesemia",  1.49,   NA,      10,      "mg/dL",        NA,
  "Hypomagnesemia",  1.5,    NA,      10,      "mg/dL",        NA,
  # Unit missing cannot grade
  "Hypomagnesemia",  1.5,    1.5,     10,           NA,        NA,
  # AVAL missing cannot grade
  "Hypomagnesemia",  NA,     1.5,     10,      "mg/dL",        NA,
)
input_magnd_cv <- expected_magnd_cv %>%
  select(-ATOXGRL)

## Test 77a: CTCAEv4 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 77a: CTCAEv4 Hypomagnesemia (SI unit)", {
  actual_magnd_si <- derive_var_atoxgr_dir(
    input_magnd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_si,
    compare = actual_magnd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 77b: CTCAEv4 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 77b: CTCAEv4 Hypomagnesemia (USCV unit)", {
  actual_magnd_cv <- derive_var_atoxgr_dir(
    input_magnd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_cv,
    compare = actual_magnd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 78a: CTCAEv5 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 78a: CTCAEv5 Hypomagnesemia (SI unit)", {
  actual_magnd_si <- derive_var_atoxgr_dir(
    input_magnd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_si,
    compare = actual_magnd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 78b: CTCAEv5 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 78b: CTCAEv5 Hypomagnesemia (USCV unit)", {
  actual_magnd_cv <- derive_var_atoxgr_dir(
    input_magnd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_cv,
    compare = actual_magnd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 78c: CTCAEv6 Hypomagnesemia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 78c: CTCAEv6 Hypomagnesemia (SI unit)", {
  actual_magnd_si <- derive_var_atoxgr_dir(
    input_magnd_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_si,
    compare = actual_magnd_si,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 78d: CTCAEv6 Hypomagnesemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 78d: CTCAEv6 Hypomagnesemia (USCV unit)", {
  actual_magnd_cv <- derive_var_atoxgr_dir(
    input_magnd_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd_cv,
    compare = actual_magnd_cv,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hyponatremia (NCICTCAEv4)
### NCICTCAEv4 and NCICTCAEv5 essentially the same (slightly different text)
### SI unit and CV unit is mmol/L
### Grade 4: <120 mmol/L
### Grade 3: <130 - 120 mmol/L
### Grade 1: <LLN - 130 mmol/L

expected_natrd <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,   ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",    119,    140,      100,    "mmol/L",  NA,
  NA_character_,   119,    140,      100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hyponatremia",  119,    140,      100,    "mmol/L",  "4",
  "Hyponatremia",  120,    140,      100,    "mmol/L",  "3",
  "Hyponatremia",  129,    140,      100,    "mmol/L",  "3",
  "Hyponatremia",  130,    140,      100,    "mmol/L",  "1",
  "Hyponatremia",  139,    140,      100,    "mmol/L",  "1",
  "Hyponatremia",  140,    140,      100,    "mmol/L",  "0",
  # ANRLO missing - can grade 3-4
  "Hyponatremia",  119,    NA,       100,    "mmol/L",  "4",
  "Hyponatremia",  120,    NA,       100,    "mmol/L",  "3",
  "Hyponatremia",  129,    NA,       100,    "mmol/L",  "3",
  # ANRLO missing - can NOT grade 0 or 1
  "Hyponatremia",  130,    NA,       100,    "mmol/L",  NA,
  "Hyponatremia",  139,    NA,       100,    "mmol/L",  NA,
  "Hyponatremia",  140,    NA,       100,    "mmol/L",  NA,
  # Unit missing cannot grade
  "Hyponatremia",  140,    140,      100,    NA,        NA,
  # AVAL missing cannot grade
  "Hyponatremia",  NA,     140,      100,    "mmol/L",  NA,
)
input_natrd <- expected_natrd %>%
  select(-ATOXGRL)

## Test 79a: CTCAEv4 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 79a: CTCAEv4 Hyponatremia (SI unit)", {
  actual_natrd <- derive_var_atoxgr_dir(
    input_natrd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natrd,
    compare = actual_natrd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 79b: CTCAEv4 Hyponatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 79b: CTCAEv4 Hyponatremia (USCV unit)", {
  actual_natrd <- derive_var_atoxgr_dir(
    input_natrd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natrd,
    compare = actual_natrd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 80a: CTCAEv5 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 80a: CTCAEv5 Hyponatremia (SI unit)", {
  actual_natrd <- derive_var_atoxgr_dir(
    input_natrd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natrd,
    compare = actual_natrd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 80b: CTCAEv5 Hyponatremia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 80b: CTCAEv5 Hyponatremia (USCV unit)", {
  actual_natrd <- derive_var_atoxgr_dir(
    input_natrd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natrd,
    compare = actual_natrd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 80c: CTCAEv6 Hyponatremia (SI unit) ----
test_that("derive_var_atoxgr_dir Test 80c: CTCAEv6 Hyponatremia (SI unit)", {
  actual_natrd <- derive_var_atoxgr_dir(
    input_natrd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_natrd,
    compare = actual_natrd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
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
  expected_phosd <- tibble::tribble(
    ~ATOXDSCL,           ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",        0.29,   1,      100,    "mmol/L",  NA,
    NA_character_,       0.29,   1,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypophosphatemia",  0.29,   1,      100,    "mmol/L",  "4",
    "Hypophosphatemia",  0.3,    1,      100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.59,   1,      100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.6,    1,      100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.79,   1,      100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.8,    1,      100,    "mmol/L",  "1",
    "Hypophosphatemia",  0.9,    1,      100,    "mmol/L",  "1",
    "Hypophosphatemia",  1,      1,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 3-4
    "Hypophosphatemia",  0.29,   NA,     100,    "mmol/L",  "4",
    "Hypophosphatemia",  0.3,    NA,     100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.59,   NA,     100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.6,    NA,     100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.79,   NA,     100,    "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypophosphatemia",  0.8,    NA,     100,    "mmol/L",  NA,
    "Hypophosphatemia",  0.9,    NA,     100,    "mmol/L",  NA,
    "Hypophosphatemia",  1,      NA,     100,    "mmol/L",  NA,
    # Unit missing cannot grade
    "Hypophosphatemia",  1,      1,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypophosphatemia",  NA,     1,      100,    "mmol/L",  NA,
  )
  input_phosd <- expected_phosd %>%
    select(-ATOXGRL)

  actual_phosd <- derive_var_atoxgr_dir(
    input_phosd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_phosd,
    compare = actual_phosd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Test 81b: CTCAEv4 Hypophosphatemia (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 81: CTCAEv4 Hypophosphatemia (USCV unit)", {
  expected_phosd <- tibble::tribble(
    ~ATOXDSCL,           ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL,
    "Not a term",        0.99,   3,       100,     "mg/dL",        NA,
    NA_character_,       0.99,   3,       100,     "mg/dL",        NA,
    # ANRLO not missing
    "Hypophosphatemia",  0.99,   3,       100,     "mg/dL",       "4",
    "Hypophosphatemia",  1,      3,       100,     "mg/dL",       "3",
    "Hypophosphatemia",  1.9,    3,       100,     "mg/dL",       "3",
    "Hypophosphatemia",  2,      3,       100,     "mg/dL",       "2",
    "Hypophosphatemia",  2.49,   3,       100,     "mg/dL",       "2",
    "Hypophosphatemia",  2.5,    3,       100,     "mg/dL",       "1",
    "Hypophosphatemia",  2.9,    3,       100,     "mg/dL",       "1",
    "Hypophosphatemia",  3,      3,       100,     "mg/dL",       "0",
    # ANRLO missing - can grade 3-4
    "Hypophosphatemia",  0.99,   NA,      100,     "mg/dL",       "4",
    "Hypophosphatemia",  1,      NA,      100,     "mg/dL",       "3",
    "Hypophosphatemia",  1.9,    NA,      100,     "mg/dL",       "3",
    "Hypophosphatemia",  2,      NA,      100,     "mg/dL",       "2",
    "Hypophosphatemia",  2.49,   NA,      100,     "mg/dL",       "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypophosphatemia",  2.5,    NA,      100,     "mg/dL",        NA,
    "Hypophosphatemia",  2.9,    NA,      100,     "mg/dL",        NA,
    "Hypophosphatemia",  3,      NA,      100,     "mg/dL",        NA,
    # Unit missing cannot grade
    "Hypophosphatemia",  3,      1,       100,          NA,        NA,
    # AVAL missing cannot grade
    "Hypophosphatemia",  NA,     1,       100,     "mg/dL",        NA,
  )
  input_phosd <- expected_phosd %>%
    select(-ATOXGRL)

  actual_phosd <- derive_var_atoxgr_dir(
    input_phosd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_phosd,
    compare = actual_phosd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# DAIDS grading----

### Acidosis
### Grade 4: pH < 7.3 with lifethreatening consequences
### Grade 3: pH < 7.3 without lifethreatening consequences
### Grade 2: pH >= 7.3 to < LLN

## Test 82: DAIDS Acidosis ----
test_that("derive_var_atoxgr_dir Test 82: DAIDS Acidosis", {
  expected_acido_daids <- tibble::tribble(
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

  input_acido_daids <- expected_acido_daids %>%
    select(-ATOXGRL)

  actual_acido_daids <- derive_var_atoxgr_dir(
    input_acido_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_acido_daids,
    compare = actual_acido_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

### Albumin, Low
### SI unit is g/L CV unit is g/dL
### Grade 3: < 20
### Grade 2: >= 20 to < 30
### Grade 1: 30 to < LLN

expected_albl_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",   35,     40,     "g/L",  NA,       1,
  NA_character_,  35,     40,     "g/L",  NA,       2,
  # ANRLO not missing
  "Albumin, Low", 19,     40,     "g/L",  "3",      3,
  "Albumin, Low", 20,     40,     "g/L",  "2",      4,
  "Albumin, Low", 29,     40,     "g/L",  "2",      5,
  "Albumin, Low", 30,     40,     "g/L",  "1",      6,
  "Albumin, Low", 39,     40,     "g/L",  "1",      7,
  "Albumin, Low", 40,     40,     "g/L",  "0",      8,
  # ANRLO missing - can grade 2 and 3
  "Albumin, Low", 19,     NA,     "g/L",  "3",      9,
  "Albumin, Low", 20,     NA,     "g/L",  "2",      10,
  "Albumin, Low", 29,     NA,     "g/L",  "2",      11,
  # ANRLO missing - can NOT grade 0 or 1
  "Albumin, Low", 30,     NA,     "g/L",  NA,       12,
  "Albumin, Low", 39,     NA,     "g/L",  NA,       13,
  "Albumin, Low", 40,     NA,     "g/L",  NA,       14,
  # AVALU missing cannot grade
  "Albumin, Low", 40,     40,     NA,     NA,       15,
  # AVAL missing cannot grade
  "Albumin, Low", NA,     40,     "g/L",  NA,       16,
)

input_albl_daids_si <- expected_albl_daids_si %>%
  select(-ATOXGRL)

expected_albl_daids_cv <- expected_albl_daids_si %>%
  mutate(
    AVAL = AVAL / 10,
    ANRLO = ANRLO / 10,
    AVALU = if_else(is.na(AVALU), NA_character_, "g/dL")
  )

input_albl_daids_cv <- expected_albl_daids_cv %>%
  select(-ATOXGRL)

## Test 83a: DAIDS Albumin, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 83a: DAIDS Albumin, Low (SI unit)", {
  actual_albl_daids_si <- derive_var_atoxgr_dir(
    input_albl_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albl_daids_si,
    compare = actual_albl_daids_si,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 83b: DAIDS Albumin, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 83b: DAIDS Albumin, Low (USCV unit)", {
  actual_albl_daids_cv <- derive_var_atoxgr_dir(
    input_albl_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albl_daids_cv,
    compare = actual_albl_daids_cv,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})


### Alkaline Phosphatase, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 84: DAIDS Alkaline Phosphatase, High ----
test_that("derive_var_atoxgr_dir Test 84: DAIDS Alkaline Phosphatase, High", {
  expected_alkphi_daids <- tibble::tribble(
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

  input_alkphi_daids <- expected_alkphi_daids %>%
    select(-ATOXGRH)

  actual_alkphi_daids <- derive_var_atoxgr_dir(
    input_alkphi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkphi_daids,
    compare = actual_alkphi_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### Alkalosis
### Grade 4: pH > 7.5 with lifethreatening consequences
### Grade 3: pH > 7.5 without lifethreatening consequences
### Grade 2: pH > ULN to <= 7.5

## Test 85: DAIDS Alkalosis ----
test_that("derive_var_atoxgr_dir Test 85: DAIDS Alkalosis", {
  expected_alkalo_daids <- tibble::tribble(
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

  input_alkalo_daids <- expected_alkalo_daids %>%
    select(-ATOXGRH)

  actual_alkalo_daids <- derive_var_atoxgr_dir(
    input_alkalo_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkalo_daids,
    compare = actual_alkalo_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})


### ALT, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 86: DAIDS ALT, High ----
test_that("derive_var_atoxgr_dir Test 86: DAIDS ALT, High", {
  expected_alti_daids <- tibble::tribble(
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

  input_alti_daids <- expected_alti_daids %>%
    select(-ATOXGRH)

  actual_alti_daids <- derive_var_atoxgr_dir(
    input_alti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alti_daids,
    compare = actual_alti_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})


### Amylase, High
### Grade 4: >= 5.0 x ULN
### Grade 3: 3.0 to < 5.0 x ULN
### Grade 2: 1.5 to < 3.0 x ULN
### Grade 1: 1.1 to < 1.5 x ULN

## Test 87: DAIDS Amylase, High ----
test_that("derive_var_atoxgr_dir Test 87: DAIDS Amylase, High", {
  expected_amyli_daids <- tibble::tribble(
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

  input_amyli_daids <- expected_amyli_daids %>%
    select(-ATOXGRH)

  actual_amyli_daids <- derive_var_atoxgr_dir(
    input_amyli_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_amyli_daids,
    compare = actual_amyli_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### AST, High
### Grade 4: >= 10.0 x ULN
### Grade 3: 5.0 to < 10.0 x ULN
### Grade 2: 2.5 to < 5.0 x ULN
### Grade 1: 1.25 to < 2.5 x ULN

## Test 88: DAIDS AST, High ----
test_that("derive_var_atoxgr_dir Test 88: DAIDS AST, High", {
  expected_asti_daids <- tibble::tribble(
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

  input_asti_daids <- expected_asti_daids %>%
    select(-ATOXGRH)

  actual_asti_daids <- derive_var_atoxgr_dir(
    input_asti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_asti_daids,
    compare = actual_asti_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})


### Bicarbonate, Low
### SI + cv unit is mmol/L
### Grade 4: < 8.0 mmol/L
### Grade 3: 8.0 -< 11.0 mmol/L
### Grade 2: 11.0 -< 16.0 mmol/L
### Grade 1: 16.0 mmol/L -< LLN

expected_bicarbd_daids <- tibble::tribble(
  ~ATOXDSCL,           ~AVAL,  ~ANRLO, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",        22,     20,     "mmol/L",  NA,       1,
  NA_character_,       22,     20,     "mmol/L",  NA,       2,
  # ANRLO not missing
  "Bicarbonate, Low",  7.9,    20,     "mmol/L",  "4",      3,
  "Bicarbonate, Low",  8,      20,     "mmol/L",  "3",      4,
  "Bicarbonate, Low",  10.9,   20,     "mmol/L",  "3",      5,
  "Bicarbonate, Low",  11,     20,     "mmol/L",  "2",      6,
  "Bicarbonate, Low",  15.9,   20,     "mmol/L",  "2",      7,
  "Bicarbonate, Low",  16,     20,     "mmol/L",  "1",      8,
  "Bicarbonate, Low",  19,     20,     "mmol/L",  "1",      9,
  "Bicarbonate, Low",  20,     20,     "mmol/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Bicarbonate, Low",  7.9,    NA,     "mmol/L",  "4",      11,
  "Bicarbonate, Low",  8,      NA,     "mmol/L",  "3",      12,
  "Bicarbonate, Low",  10.9,   NA,     "mmol/L",  "3",      13,
  "Bicarbonate, Low",  11,     NA,     "mmol/L",  "2",      14,
  "Bicarbonate, Low",  15.9,   NA,     "mmol/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Bicarbonate, Low",  16,     NA,     "mmol/L",  NA,       16,
  "Bicarbonate, Low",  19,     NA,     "mmol/L",  NA,       17,
  "Bicarbonate, Low",  20,     NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Bicarbonate, Low",  20,     20,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Bicarbonate, Low",  NA,     20,     "mmol/L",  NA,       20,
)
input_bicarbd_daids <- expected_bicarbd_daids %>%
  select(-ATOXGRL)

## Test 89a: DAIDS Bicarbonate, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 89a: DAIDS Bicarbonate, Low (SI unit)", {
  actual_bicarbd_daids <- derive_var_atoxgr_dir(
    input_bicarbd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bicarbd_daids,
    compare = actual_bicarbd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 89b: DAIDS Bicarbonate, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 89b: DAIDS Bicarbonate, Low (USCV unit)", {
  actual_bicarbd_daids <- derive_var_atoxgr_dir(
    input_bicarbd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bicarbd_daids,
    compare = actual_bicarbd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})


### Direct Bilirubin, High

### SI unit is umol/L
### CV unit is mg/dL
### 17.1 used as conversion from "mg/dL" to "umol/L"

### > 28 days of age

### Grade 4: > ULN

## SI unit
expected_dbiligt28d_daids_si <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",               7,      8,      "umol/L", NA,       1,
  NA_character_,              7,      8,      "umol/L", NA,       2,
  # ANRHI not missing
  "Direct Bilirubin, High",   8.1,    8,      "umol/L", "4",      3,
  "Direct Bilirubin, High",   8,      8,      "umol/L", "0",      4,
  "Direct Bilirubin, High",   7.9,    8,      "umol/L", "0",      5,
  # ANRHI missing cannot still grade 0 or 4
  "Direct Bilirubin, High",   8.1,    NA,     "umol/L", NA,       6,
  "Direct Bilirubin, High",   8,      NA,     "umol/L", NA,       7,
  "Direct Bilirubin, High",   7.9,    NA,     "umol/L", NA,       8,
  # AVAL missing cannot grade
  "Direct Bilirubin, High",   NA,     8,      "umol/L", NA,       9,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-30")
  )

## CV unit
expected_dbiligt28d_daids_cv <- expected_dbiligt28d_daids_si %>%
  mutate(
    AVALU = "mg/dL",
    TESTNUM = TESTNUM + 30
  )

### <= 28 days of age

### Grade 4: > 2 mg/dL (> 34.2 umol/L)
### Grade 3: > 1.5 to <= 2 mg/dL (> 25.65 to <= 34.2 umol/L)
### Grade 2: > 1 to <= 1.5 mg/dL (> 17.1 to <= 25.65 umol/L)
### Grade 1: ULN to <= 1 mg/dL (ULN to <= 17.1 umol/L)

## SI unit
expected_dbilile28d_daids_si <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",               7,      8,      "umol/L", NA,       10,
  NA_character_,              7,      8,      "umol/L", NA,       11,
  # ANRHI not missing
  "Direct Bilirubin, High",   34.3,   8,      "umol/L", "4",      12,
  "Direct Bilirubin, High",   34.2,   8,      "umol/L", "3",      13,
  "Direct Bilirubin, High",   25.66,  8,      "umol/L", "3",      14,
  "Direct Bilirubin, High",   25.65,  8,      "umol/L", "2",      15,
  "Direct Bilirubin, High",   17.19,  8,      "umol/L", "2",      16,
  "Direct Bilirubin, High",   17.1,   8,      "umol/L", "1",      17,
  "Direct Bilirubin, High",   8,      8,      "umol/L", "1",      18,
  "Direct Bilirubin, High",   7.9,    8,      "umol/L", "0",      19,
  # ANRHI missing can still grade 2 - 4
  "Direct Bilirubin, High",   34.3,   NA,     "umol/L", "4",      20,
  "Direct Bilirubin, High",   34.2,   NA,     "umol/L", "3",      21,
  "Direct Bilirubin, High",   25.66,  NA,     "umol/L", "3",      22,
  "Direct Bilirubin, High",   25.65,  NA,     "umol/L", "2",      23,
  "Direct Bilirubin, High",   17.19,  NA,     "umol/L", "2",      24,
  # ANRHI missing cannot still grade 0 - 1
  "Direct Bilirubin, High",   17.1,   NA,     "umol/L", NA,       25,
  "Direct Bilirubin, High",   8,      NA,     "umol/L", NA,       26,
  "Direct Bilirubin, High",   7.9,    NA,     "umol/L", NA,       27,
  # AVAL missing cannot grade
  "Direct Bilirubin, High",   NA,     8,      "umol/L", NA,       28,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-29")
  )

## CV unit
expected_dbilile28d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,  ~ATOXGRH, ~TESTNUM,
  "Not a term",               2.1,    0.5,    "mg/dL", NA,       10,
  NA_character_,              2.1,    0.5,    "mg/dL", NA,       11,
  # ANRHI not missing
  "Direct Bilirubin, High",   2.1,    0.5,    "mg/dL", "4",      12,
  "Direct Bilirubin, High",   2,      0.5,    "mg/dL", "3",      13,
  "Direct Bilirubin, High",   1.51,   0.5,    "mg/dL", "3",      14,
  "Direct Bilirubin, High",   1.5,    0.5,    "mg/dL", "2",      15,
  "Direct Bilirubin, High",   1.1,    0.5,    "mg/dL", "2",      16,
  "Direct Bilirubin, High",   1,      0.5,    "mg/dL", "1",      17,
  "Direct Bilirubin, High",   0.5,    0.5,    "mg/dL", "1",      18,
  "Direct Bilirubin, High",   0.49,   0.5,    "mg/dL", "0",      19,
  # ANRHI missing can still grade 2 - 4
  "Direct Bilirubin, High",   2.1,    NA,     "mg/dL", "4",      20,
  "Direct Bilirubin, High",   2,      NA,     "mg/dL", "3",      21,
  "Direct Bilirubin, High",   1.51,   NA,     "mg/dL", "3",      22,
  "Direct Bilirubin, High",   1.5,    NA,     "mg/dL", "2",      23,
  "Direct Bilirubin, High",   1.1,    NA,     "mg/dL", "2",      24,
  # ANRHI missing cannot still grade 0 - 1
  "Direct Bilirubin, High",   1,      NA,     "mg/dL", NA,       25,
  "Direct Bilirubin, High",   0.5,    NA,     "mg/dL", NA,       26,
  "Direct Bilirubin, High",   0.49,   NA,     "mg/dL", NA,       27,
  # AVAL missing cannot grademg/dL
  "Direct Bilirubin, High",   NA,     0.5,    "mg/dL", NA,       28,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-29")
  )


### add subjects with missing ADT or BRTHDT

## SI unit
expected_dbilinoage_daids_si <- expected_dbilile28d_daids_si %>%
  filter(TESTNUM %in% c(18, 19)) %>%
  mutate(
    ADT = if_else(TESTNUM == 18, NA, ADT),
    BRTHDT = if_else(TESTNUM == 19, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 18, 29, 30)
  )

## CV unit
expected_dbilinoage_daids_cv <- expected_dbilile28d_daids_cv %>%
  filter(TESTNUM %in% c(48, 49)) %>%
  mutate(
    ADT = if_else(TESTNUM == 48, NA, ADT),
    BRTHDT = if_else(TESTNUM == 49, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 48, 59, 60)
  )

### put all SI data together
expected_dbili_daids_si <- expected_dbilinoage_daids_si %>%
  bind_rows(
    expected_dbilile28d_daids_si,
    expected_dbiligt28d_daids_si
  )

### put all CV data together
expected_dbili_daids_cv <- expected_dbilinoage_daids_cv %>%
  bind_rows(
    expected_dbilile28d_daids_cv,
    expected_dbiligt28d_daids_cv
  )

## Test 90a: DAIDS Direct Bilirubin, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 90: DAIDS Direct Bilirubin, High (SI unit)", {
  input_dbili_daids_si <- expected_dbili_daids_si %>%
    select(-ATOXGRH)

  actual_dbili_daids_si <- derive_var_atoxgr_dir(
    input_dbili_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_dbili_daids_si,
    compare = actual_dbili_daids_si,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 90b: DAIDS Direct Bilirubin, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 90b: DAIDS Direct Bilirubin, High (USCV unit)", {
  input_dbili_daids_cv <- expected_dbili_daids_cv %>%
    select(-ATOXGRH)

  actual_dbili_daids_cv <- derive_var_atoxgr_dir(
    input_dbili_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_dbili_daids_cv,
    compare = actual_dbili_daids_cv,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### Total Bilirubin, High

### > 28 days of age

### Grade 4: >= 5.0 x ULN
### Grade 3: 2.6 to < 5.0 x ULN
### Grade 2: 1.6 to < 2.6 x ULN
### Grade 1: 1.1 to < 1.6 x ULN

expected_tbiligt28d_daids <- tibble::tribble(
  ~ATOXDSCH,                ~AVAL,  ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",             9,      10,     "umol/L",  NA,       1,
  NA_character_,            9,      10,     "umol/L",  NA,       2,
  # ANRHI not missing
  "Total Bilirubin, High",  50,     10,     "umol/L",  "4",      3,
  "Total Bilirubin, High",  49,     10,     "umol/L",  "3",      4,
  "Total Bilirubin, High",  26,     10,     "umol/L",  "3",      5,
  "Total Bilirubin, High",  25.9,   10,     "umol/L",  "2",      6,
  "Total Bilirubin, High",  16,     10,     "umol/L",  "2",      7,
  "Total Bilirubin, High",  15.9,   10,     "umol/L",  "1",      8,
  "Total Bilirubin, High",  11,     10,     "umol/L",  "1",      9,
  "Total Bilirubin, High",  10.9,   10,     "umol/L",  "0",      10,
  # Unit missing can grade - grade based on comparison of AVAL with ANRHI
  "Total Bilirubin, High",  10.9,   10,     NA,        "0",      11,
  # ANRHI missing - cannot grade
  "Total Bilirubin, High",  10.9,   NA,     "umol/L",  NA,       12,
  # AVAL missing cannot grade
  "Total Bilirubin, High",  NA,     10,     "umol/L",  NA,       13,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-30")
  )

### make Age <= 28 all results NA for ATOXGRH
expected_tbilile28d_daids <- expected_tbiligt28d_daids %>%
  mutate(
    ADT = lubridate::ymd("2023-01-29"),
    ATOXGRH = NA_character_,
    TESTNUM = TESTNUM + 13
  )

### make Age missing results NA for ATOXGRH
expected_tbilinoage_daids <- expected_tbiligt28d_daids %>%
  filter(TESTNUM %in% c(10, 11)) %>%
  mutate(
    ADT = if_else(TESTNUM == 10, NA, ADT),
    BRTHDT = if_else(TESTNUM == 11, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 10, 27, 28)
  )

expected_tbili_daids <- expected_tbilinoage_daids %>%
  bind_rows(
    expected_tbiligt28d_daids,
    expected_tbilile28d_daids
  )


## Test 91: DAIDS Total Bilirubin, High ----
test_that("derive_var_atoxgr_dir Test 91: DAIDS Total Bilirubin, High", {
  input_tbili_daids <- expected_tbili_daids %>%
    select(-ATOXGRH)

  actual_tbili_daids <- derive_var_atoxgr_dir(
    input_tbili_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_tbili_daids,
    compare = actual_tbili_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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

expected_calcige7d_daids_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    "mmol/L",  NA,       1,
  NA_character_,   3.5,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Calcium, High", 3.38,   "mmol/L",  "4",      3,
  "Calcium, High", 3.37,   "mmol/L",  "3",      4,
  "Calcium, High", 3.13,   "mmol/L",  "3",      5,
  "Calcium, High", 3.12,   "mmol/L",  "2",      6,
  "Calcium, High", 2.88,   "mmol/L",  "2",      7,
  "Calcium, High", 2.87,   "mmol/L",  "1",      8,
  "Calcium, High", 2.65,   "mmol/L",  "1",      9,
  "Calcium, High", 2.64,   "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Calcium, High", 2.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mmol/L",  NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

expected_calcige7d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",     13.5,   "mg/dL",        NA,  1,
  NA_character_,    13.5,   "mg/dL",        NA,  2,
  # ANRHI not missing
  "Calcium, High",  13.5,   "mg/dL",       "4",  3,
  "Calcium, High",  13.49,  "mg/dL",       "3",  4,
  "Calcium, High",  12.5,   "mg/dL",       "3",  5,
  "Calcium, High",  12.49,  "mg/dL",       "2",  6,
  "Calcium, High",  11.5,   "mg/dL",       "2",  7,
  "Calcium, High",  11.49,  "mg/dL",       "1",  8,
  "Calcium, High",  10.6,   "mg/dL",       "1",  9,
  "Calcium, High",  10.59,  "mg/dL",       "0",  10,
  # Unit missing cannot grade
  "Calcium, High",  10.5,        NA,        NA,  11,
  # AVAL missing cannot grade
  "Calcium, High",  NA,     "mg/dL",        NA,  12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

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

expected_calcilt7d_daids_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    "mmol/L",  NA,       13,
  NA_character_,   3.5,    "mmol/L",  NA,       14,
  # ANRHI not missing
  "Calcium, High", 3.38,   "mmol/L",  "4",      15,
  "Calcium, High", 3.37,   "mmol/L",  "3",      16,
  "Calcium, High", 3.23,   "mmol/L",  "3",      17,
  "Calcium, High", 3.22,   "mmol/L",  "2",      18,
  "Calcium, High", 3.1,    "mmol/L",  "2",      19,
  "Calcium, High", 3.09,   "mmol/L",  "1",      20,
  "Calcium, High", 2.88,   "mmol/L",  "1",      21,
  "Calcium, High", 2.87,   "mmol/L",  "0",      22,
  # Unit missing cannot grade
  "Calcium, High", 3.5,    NA,        NA,       23,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mmol/L",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

expected_calcilt7d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",    13.5,   "mg/dL",  NA,       13,
  NA_character_,   13.5,   "mg/dL",  NA,       14,
  # ANRHI not missing
  "Calcium, High", 13.5,   "mg/dL",  "4",      15,
  "Calcium, High", 13.49,  "mg/dL",  "3",      16,
  "Calcium, High", 12.9,   "mg/dL",  "3",      17,
  "Calcium, High", 12.89,  "mg/dL",  "2",      18,
  "Calcium, High", 12.4,   "mg/dL",  "2",      19,
  "Calcium, High", 12.39,  "mg/dL",  "1",      20,
  "Calcium, High", 11.5,   "mg/dL",  "1",      21,
  "Calcium, High", 11.49,  "mg/dL",  "0",      22,
  # Unit missing cannot grade
  "Calcium, High", 11.49,  NA,       NA,       23,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mg/dL",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

## create data with age missing for SI unit
expected_calcinoage_daids_si <- expected_calcige7d_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## create data with age missing for CV unit
expected_calcinoage_daids_cv <- expected_calcige7d_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

# put all SI data together
expected_calci_daids_si <- expected_calcinoage_daids_si %>%
  bind_rows(
    expected_calcige7d_daids_si,
    expected_calcilt7d_daids_si
  )

# put all CV data together
expected_calci_daids_cv <- expected_calcinoage_daids_cv %>%
  bind_rows(
    expected_calcige7d_daids_cv,
    expected_calcilt7d_daids_cv
  )
## Test 92a: DAIDS Calcium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 92a: DAIDS Calcium, High (SI unit)", {
  input_calci_daids_si <- expected_calci_daids_si %>%
    select(-ATOXGRH)

  actual_calci_daids_si <- derive_var_atoxgr_dir(
    input_calci_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_daids_si,
    compare = actual_calci_daids_si,
    keys = c("TESTNUM")
  )
})

## Test 92b: DAIDS Calcium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 92b: DAIDS Calcium, High (USCV unit)", {
  input_calci_daids_cv <- expected_calci_daids_cv %>%
    select(-ATOXGRH)

  actual_calci_daids_cv <- derive_var_atoxgr_dir(
    input_calci_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci_daids_cv,
    compare = actual_calci_daids_cv,
    keys = c("TESTNUM")
  )
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

expected_calioni_daids_si <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",              1.8,   1.1,    1.4,    "mmol/L", NA,       1,
  NA_character_,             1.79,  1.1,    1.4,    "mmol/L", NA,       2,
  # ANRHI not missing
  "Calcium (Ionized), High", 1.8,   1.1,    1.4,    "mmol/L", "4",      3,
  "Calcium (Ionized), High", 1.79,  1.1,    1.4,    "mmol/L", "3",      4,
  "Calcium (Ionized), High", 1.6,   1.1,    1.4,    "mmol/L", "3",      5,
  "Calcium (Ionized), High", 1.59,  1.1,    1.4,    "mmol/L", "2",      6,
  "Calcium (Ionized), High", 1.5,   1.1,    1.4,    "mmol/L", "2",      7,
  "Calcium (Ionized), High", 1.49,  1.1,    1.4,    "mmol/L", "1",      8,
  "Calcium (Ionized), High", 1.41,  1.1,    1.4,    "mmol/L", "1",      9,
  "Calcium (Ionized), High", 1.4,   1.1,    1.4,    "mmol/L", "0",      10,
  # ANRHI missing - can grade 2-4
  "Calcium (Ionized), High", 1.8,   1.1,    NA,     "mmol/L", "4",      11,
  "Calcium (Ionized), High", 1.79,  1.1,    NA,     "mmol/L", "3",      12,
  "Calcium (Ionized), High", 1.6,   1.1,    NA,     "mmol/L", "3",      13,
  "Calcium (Ionized), High", 1.59,  1.1,    NA,     "mmol/L", "2",      14,
  "Calcium (Ionized), High", 1.5,   1.1,    NA,     "mmol/L", "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Calcium (Ionized), High", 1.49,  1.1,    NA,     "mmol/L", NA,       16,
  "Calcium (Ionized), High", 1.41,  1.1,    NA,     "mmol/L", NA,       17,
  "Calcium (Ionized), High", 1.4,   1.1,    NA,     "mmol/L", NA,       18,
  # Unit missing cannot grade
  "Calcium (Ionized), High", 1.3,   1.1,    1.4,    NA,       NA,       19,
  # AVAL missing cannot grade
  "Calcium (Ionized), High", NA,    1.1,    1.4,    "mmol/L", NA,       20,
)
input_calioni_daids_si <- expected_calioni_daids_si %>%
  select(-ATOXGRH)

expected_calioni_daids_cv <- expected_calioni_daids_si %>%
  mutate(
    AVAL = 4 * AVAL,
    ANRLO = 4 * ANRLO,
    ANRHI = 4 * ANRHI,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )
input_calioni_daids_cv <- expected_calioni_daids_cv %>%
  select(-ATOXGRH)

## Test 93a: DAIDS Calcium (Ionized), High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 93a: DAIDS Calcium (Ionized), High (SI unit)", {
  actual_calioni_daids_si <- derive_var_atoxgr_dir(
    input_calioni_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_daids_si,
    compare = actual_calioni_daids_si,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 93b: DAIDS Calcium (Ionized), High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 93b: DAIDS Calcium (Ionized), High (USCV unit)", {
  actual_calioni_daids_cv <- derive_var_atoxgr_dir(
    input_calioni_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni_daids_cv,
    compare = actual_calioni_daids_cv,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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

expected_calcdge7d_daids_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",    2.2,    "mmol/L",  NA,       1,
  NA_character_,   2.2,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Calcium, Low",  1.52,   "mmol/L",  "4",      3,
  "Calcium, Low",  1.53,   "mmol/L",  "3",      4,
  "Calcium, Low",  1.74,   "mmol/L",  "3",      5,
  "Calcium, Low",  1.75,   "mmol/L",  "2",      6,
  "Calcium, Low",  1.94,   "mmol/L",  "2",      7,
  "Calcium, Low",  1.95,   "mmol/L",  "1",      8,
  "Calcium, Low",  2.09,   "mmol/L",  "1",      9,
  "Calcium, Low",  2.1,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Calcium, Low",  2.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mmol/L",  NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

expected_calcdge7d_daids_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  "Not a term",    6.09,   "mg/dL",        NA,  1,
  NA_character_,   6.09,   "mg/dL",        NA,  2,
  # ANRLO not missing
  "Calcium, Low",  6.09,   "mg/dL",       "4",  3,
  "Calcium, Low",  6.1,    "mg/dL",       "3",  4,
  "Calcium, Low",  6.9,    "mg/dL",       "3",  5,
  "Calcium, Low",  7,      "mg/dL",       "2",  6,
  "Calcium, Low",  7.7,    "mg/dL",       "2",  7,
  "Calcium, Low",  7.8,    "mg/dL",       "1",  8,
  "Calcium, Low",  8.3,    "mg/dL",       "1",  9,
  "Calcium, Low",  8.4,    "mg/dL",       "0",  10,
  # Unit missing cannot grade
  "Calcium, Low",  8.4,         NA,        NA,  11,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mg/dL",        NA,  12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

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
expected_calcdlt7d_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   2.2,   "mmol/L", NA,       13,
  NA_character_,  2.2,   "mmol/L", NA,       14,
  # ANRLO not missing
  "Calcium, Low", 1.37,  "mmol/L", "4",      15,
  "Calcium, Low", 1.38,  "mmol/L", "3",      16,
  "Calcium, Low", 1.49,  "mmol/L", "3",      17,
  "Calcium, Low", 1.5,   "mmol/L", "2",      18,
  "Calcium, Low", 1.62,  "mmol/L", "2",      19,
  "Calcium, Low", 1.63,  "mmol/L", "1",      20,
  "Calcium, Low", 1.87,  "mmol/L", "1",      21,
  "Calcium, Low", 1.88,  "mmol/L", "0",      22,
  # Unit missing cannot grade
  "Calcium, Low", 2.2,   NA,       NA,       23,
  # AVAL missing cannot grade
  "Calcium, Low", NA,    "mmol/L", NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )
### CV unit
expected_calcdlt7d_daids_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  "Not a term",    5.49,   "mg/dL",        NA,  13,
  NA_character_,   5.49,   "mg/dL",        NA,  14,
  # ANRLO not missing
  "Calcium, Low",  5.49,   "mg/dL",       "4",  15,
  "Calcium, Low",  5.5,    "mg/dL",       "3",  16,
  "Calcium, Low",  5.99,   "mg/dL",       "3",  17,
  "Calcium, Low",  6,      "mg/dL",       "2",  18,
  "Calcium, Low",  6.49,   "mg/dL",       "2",  19,
  "Calcium, Low",  6.5,    "mg/dL",       "1",  20,
  "Calcium, Low",  7.49,   "mg/dL",       "1",  21,
  "Calcium, Low",  7.5,    "mg/dL",       "0",  22,
  # Unit missing cannot grade
  "Calcium, Low",  7.5,         NA,        NA,  23,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mg/dL",        NA,  24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

## create SI records with age missing
expected_calcdnoage_daids_si <- expected_calcdge7d_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## create CV records with age missing
expected_calcdnoage_daids_cv <- expected_calcdge7d_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## put all SI records together
expected_calcd_daids_si <- expected_calcdnoage_daids_si %>%
  bind_rows(
    expected_calcdge7d_daids_si,
    expected_calcdlt7d_daids_si
  )

## put all CV records together
expected_calcd_daids_cv <- expected_calcdnoage_daids_cv %>%
  bind_rows(
    expected_calcdge7d_daids_cv,
    expected_calcdlt7d_daids_cv
  )


## Test 94a: DAIDS Calcium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 94a: DAIDS Calcium, Low (SI unit)", {
  input_calcd_daids_si <- expected_calcd_daids_si %>%
    select(-ATOXGRL)

  actual_calcd_daids_si <- derive_var_atoxgr_dir(
    input_calcd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_daids_si,
    compare = actual_calcd_daids_si,
    keys = c("TESTNUM")
  )
})

## Test 94b: DAIDS Calcium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 94b: DAIDS Calcium, Low (USCV unit)", {
  input_calcd_daids_cv <- expected_calcd_daids_cv %>%
    select(-ATOXGRL)

  actual_calcd_daids_cv <- derive_var_atoxgr_dir(
    input_calcd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd_daids_cv,
    compare = actual_calcd_daids_cv,
    keys = c("TESTNUM")
  )
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

expected_caliond_daids_si <- tibble::tribble(
  ~ATOXDSCL,                ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",             0.79,  1.1,    100,    "mmol/L", NA,       1,
  NA_character_,            0.79,  1.1,    100,    "mmol/L", NA,       2,
  # ANRLO not missing
  "Calcium (Ionized), Low", 0.79,  1.1,    100,    "mmol/L", "4",      3,
  "Calcium (Ionized), Low", 0.8,   1.1,    100,    "mmol/L", "3",      4,
  "Calcium (Ionized), Low", 0.89,  1.1,    100,    "mmol/L", "3",      5,
  "Calcium (Ionized), Low", 0.9,   1.1,    100,    "mmol/L", "2",      6,
  "Calcium (Ionized), Low", 0.99,  1.1,    100,    "mmol/L", "2",      7,
  "Calcium (Ionized), Low", 1,     1.1,    100,    "mmol/L", "1",      8,
  "Calcium (Ionized), Low", 1.09,  1.1,    100,    "mmol/L", "1",      9,
  "Calcium (Ionized), Low", 1.1,   1.1,    100,    "mmol/L", "0",      10,
  # ANRLO missing - can grade 2-4
  "Calcium (Ionized), Low", 0.79,  NA,     100,    "mmol/L", "4",      11,
  "Calcium (Ionized), Low", 0.8,   NA,     100,    "mmol/L", "3",      12,
  "Calcium (Ionized), Low", 0.89,  NA,     100,    "mmol/L", "3",      13,
  "Calcium (Ionized), Low", 0.9,   NA,     100,    "mmol/L", "2",      14,
  "Calcium (Ionized), Low", 0.99,  NA,     100,    "mmol/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Calcium (Ionized), Low", 1,     NA,     100,    "mmol/L", NA,       16,
  "Calcium (Ionized), Low", 1.09,  NA,     100,    "mmol/L", NA,       17,
  "Calcium (Ionized), Low", 1.1,   NA,     100,    "mmol/L", NA,       18,
  # Unit missing cannot grade
  "Calcium (Ionized), Low", 1.1,   1.1,    100,    NA,       NA,       19,
  # AVAL missing cannot grade
  "Calcium (Ionized), Low", NA,    1.1,    100,    "mmol/L", NA,       20,
)
input_caliond_daids_si <- expected_caliond_daids_si %>%
  select(-ATOXGRL)

expected_caliond_daids_cv <- expected_caliond_daids_si %>%
  mutate(
    AVAL = AVAL * 4,
    ANRLO = ANRLO * 4,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )
input_caliond_daids_cv <- expected_caliond_daids_cv %>%
  select(-ATOXGRL)

## Test 95a: DAIDS Calcium (Ionized), Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 95a: DAIDS Calcium (Ionized), Low (SI unit)", {
  actual_caliond_daids_si <- derive_var_atoxgr_dir(
    input_caliond_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_daids_si,
    compare = actual_caliond_daids_si,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 95b: DAIDS Calcium (Ionized), Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 95b: DAIDS Calcium (Ionized), Low (USCV unit)", {
  actual_caliond_daids_cv <- derive_var_atoxgr_dir(
    input_caliond_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond_daids_cv,
    compare = actual_caliond_daids_cv,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})


### Creatine Kinase, High
### Grade 4: >= 20.0 x ULN
### Grade 3: 10.0 -< 20.0 x ULN
### Grade 2: 6 -< 10.0 x ULN
### Grade 1: 3 -< 6 x ULN

## Test 96: DAIDS Creatine Kinase, High ----
test_that("derive_var_atoxgr_dir Test 96: DAIDS Creatine Kinase, High", {
  expected_cki_daids <- tibble::tribble(
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

  input_cki_daids <- expected_cki_daids %>%
    select(-ATOXGRH)

  actual_cki_daids <- derive_var_atoxgr_dir(
    input_cki_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cki_daids,
    compare = actual_cki_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})


### Creatinine, High
### Grade 4: >= 3.5 x ULN or >= 2 X BASE
### Grade 3: >1.8 -< 3.5 x ULN or 1.5 -< 2 x BASE
### Grade 2: >1.3 - 1.8 x ULN or 1.3 - < 1.5 x BASE
### Grade 1: 1.1 - 1.3 x ULN

## Test 97: DAIDS Creatinine, High ----
test_that("derive_var_atoxgr_dir Test 97: DAIDS Creatinine, High", {
  expected_creati_daids <- tibble::tribble(
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

  input_creati_daids <- expected_creati_daids %>%
    select(-ATOXGRH)

  actual_creati_daids <- derive_var_atoxgr_dir(
    input_creati_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_creati_daids,
    compare = actual_creati_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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
  expected_glucfi_daids_si <- tibble::tribble(
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

  input_glucfi_daids_si <- expected_glucfi_daids_si %>%
    select(-ATOXGRH)

  actual_glucfi_daids_si <- derive_var_atoxgr_dir(
    input_glucfi_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucfi_daids_si,
    compare = actual_glucfi_daids_si,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 98b:  DAIDS Glucose Fasting, High (USCV unit)----
test_that("derive_var_atoxgr_dir Test 98b:  DAIDS Glucose Fasting, High (USCV unit)", {
  expected_glucfi_daids_cv <- tibble::tribble(
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

  input_glucfi_daids_cv <- expected_glucfi_daids_cv %>%
    select(-ATOXGRH)

  actual_glucfi_daids_cv <- derive_var_atoxgr_dir(
    input_glucfi_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucfi_daids_cv,
    compare = actual_glucfi_daids_cv,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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
  expected_glucnfi_daids_si <- tibble::tribble(
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

  input_glucnfi_daids_si <- expected_glucnfi_daids_si %>%
    select(-ATOXGRH)

  actual_glucnfi_daids_si <- derive_var_atoxgr_dir(
    input_glucnfi_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucnfi_daids_si,
    compare = actual_glucnfi_daids_si,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 99b:  DAIDS Glucose Nonfasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 99b:  DAIDS Glucose Nonfasting, High (USCV unit)", {
  expected_glucnfi_daids_cv <- tibble::tribble(
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

  input_glucnfi_daids_cv <- expected_glucnfi_daids_cv %>%
    select(-ATOXGRH)

  actual_glucnfi_daids_cv <- derive_var_atoxgr_dir(
    input_glucnfi_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucnfi_daids_cv,
    compare = actual_glucnfi_daids_cv,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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
expected_glucdge1m_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   9.5,    "mg/L",   NA,       1,
  NA_character_,  4.1,    "mmol/L", NA,       2,
  "Glucose, Low", 1.66,   "mmol/L", "4",      3,
  "Glucose, Low", 1.67,   "mmol/L", "3",      4,
  "Glucose, Low", 2.21,   "mmol/L", "3",      5,
  "Glucose, Low", 2.22,   "mmol/L", "2",      6,
  "Glucose, Low", 3.04,   "mmol/L", "2",      7,
  "Glucose, Low", 3.05,   "mmol/L", "1",      8,
  "Glucose, Low", 3.54,   "mmol/L", "1",      9,
  "Glucose, Low", 3.55,   "mmol/L", "0",      10,
  # AVALU missing cannot grade
  "Glucose, Low", 4,      NA,       NA,       11,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mmol/L", NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-30"),
  )

### CV unit
expected_glucdge1m_daids_cv <- tibble::tribble(
  ~ATOXDSCL, ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term", 29, "mg/L", NA, 1,
  NA_character_, 29, "mg/dL", NA, 2,
  "Glucose, Low", 29, "mg/dL", "4", 3,
  "Glucose, Low", 30, "mg/dL", "3", 4,
  "Glucose, Low", 39, "mg/dL", "3", 5,
  "Glucose, Low", 40, "mg/dL", "2", 6,
  "Glucose, Low", 54, "mg/dL", "2", 7,
  "Glucose, Low", 55, "mg/dL", "1", 8,
  "Glucose, Low", 64, "mg/dL", "1", 9,
  "Glucose, Low", 65, "mg/dL", "0", 10,
  # AVALU missing cannot grade
  "Glucose, Low", 65, NA, NA, 11,
  # AVAL missing cannot grade
  "Glucose, Low", NA, "mg/dL", NA, 12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-30"),
  )

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

expected_glucdlt1m_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   9.5,    "mg/L",   NA,       13,
  NA_character_,  4.1,    "mmol/L", NA,       14,
  "Glucose, Low", 1.66,   "mmol/L", "4",      15,
  "Glucose, Low", 1.67,   "mmol/L", "3",      16,
  "Glucose, Low", 2.21,   "mmol/L", "3",      17,
  "Glucose, Low", 2.22,   "mmol/L", "2",      18,
  "Glucose, Low", 2.77,   "mmol/L", "2",      19,
  "Glucose, Low", 2.78,   "mmol/L", "1",      20,
  "Glucose, Low", 2.99,   "mmol/L", "1",      21,
  "Glucose, Low", 3,      "mmol/L", "0",      22,
  # AVALU missing cannot grade
  "Glucose, Low", 4,      NA,       NA,       23,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mmol/L", NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-29"),
  )

expected_glucdlt1m_daids_cv <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   29,     "mg/L",   NA,       13,
  NA_character_,  29,     "mg/dL",  NA,       14,
  "Glucose, Low", 29,     "mg/dL",  "4",      15,
  "Glucose, Low", 30,     "mg/dL",  "3",      16,
  "Glucose, Low", 39,     "mg/dL",  "3",      17,
  "Glucose, Low", 40,     "mg/dL",  "2",      18,
  "Glucose, Low", 49,     "mg/dL",  "2",      19,
  "Glucose, Low", 50,     "mg/dL",  "1",      20,
  "Glucose, Low", 54,     "mg/dL",  "1",      21,
  "Glucose, Low", 55,     "mg/dL",  "0",      22,
  # AVALU missing cannot grade
  "Glucose, Low", 55,     NA,       NA,       23,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mg/dL",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-29"),
  )

## create SI records with age missing
expected_glucdnoage_daids_si <- expected_glucdge1m_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ATOXGRL = NA_character_,
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )
## create CV records with age missing
expected_glucdnoage_daids_cv <- expected_glucdge1m_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ATOXGRL = NA_character_,
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## put all SI records together
expected_glucd_daids_si <- expected_glucdnoage_daids_si %>%
  bind_rows(
    expected_glucdge1m_daids_si,
    expected_glucdlt1m_daids_si
  )

## put all CV records together
expected_glucd_daids_cv <- expected_glucdnoage_daids_cv %>%
  bind_rows(
    expected_glucdge1m_daids_cv,
    expected_glucdlt1m_daids_cv
  )


## Test 100a:  DAIDS Glucose, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 100a:  DAIDS Glucose, Low (SI unit)", {
  input_glucd_daids_si <- expected_glucd_daids_si %>%
    select(-ATOXGRL)

  actual_glucd_daids_si <- derive_var_atoxgr_dir(
    input_glucd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucd_daids_si,
    compare = actual_glucd_daids_si,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 100b:  DAIDS Glucose, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 100b:  DAIDS Glucose, Low (USCV unit)", {
  input_glucd_daids_cv <- expected_glucd_daids_cv %>%
    select(-ATOXGRL)

  actual_glucd_daids_cv <- derive_var_atoxgr_dir(
    input_glucd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glucd_daids_cv,
    compare = actual_glucd_daids_cv,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})


### Lactate, High

### Grade 2: >= 2.0 x ULN
### Grade 1: ULN -< 2.0 x ULN

## Test 101:  DAIDS Lactate, High ----
test_that("derive_var_atoxgr_dir Test 101:  DAIDS Lactate, High", {
  expected_lacti_daids <- tibble::tribble(
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

  input_lacti_daids <- expected_lacti_daids %>%
    select(-ATOXGRH)

  actual_lacti_daids <- derive_var_atoxgr_dir(
    input_lacti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lacti_daids,
    compare = actual_lacti_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})


### Lipase, High

### Grade 4: >= 5.0 x ULN
### Grade 3: 3.0 -< 5.0 x ULN
### Grade 2: 1.5 -< 3.0 x ULN
### Grade 1: 1.1 -< 1.5 x ULN

## Test 102:  DAIDS Lipase, High ----
test_that("derive_var_atoxgr_dir Test 102:  DAIDS Lipase, High", {
  expected_lipi_daids <- tibble::tribble(
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

  input_lipi_daids <- expected_lipi_daids %>%
    select(-ATOXGRH)

  actual_lipi_daids <- derive_var_atoxgr_dir(
    input_lipi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lipi_daids,
    compare = actual_lipi_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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
expected_cholfige18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 3.5,    "mmol/L",  NA,       1,
  NA_character_,                3.5,    "mmol/L",  NA,       2,
  "Cholesterol, Fasting, High", 7.77,   "mmol/L",  "3",      3,
  "Cholesterol, Fasting, High", 7.76,   "mmol/L",  "2",      4,
  "Cholesterol, Fasting, High", 6.19,   "mmol/L",  "2",      5,
  "Cholesterol, Fasting, High", 6.18,   "mmol/L",  "1",      6,
  "Cholesterol, Fasting, High", 5.18,   "mmol/L",  "1",      7,
  "Cholesterol, Fasting, High", 5.17,   "mmol/L",  "0",      8,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 3.5,    NA,        NA,       9,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mmol/L",  NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

### CV data
expected_cholfige18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 300,    "mg/dL",   NA,       1,
  NA_character_,                300,    "mg/dL",   NA,       2,
  "Cholesterol, Fasting, High", 300,    "mg/dL",   "3",      3,
  "Cholesterol, Fasting, High", 299,    "mg/dL",   "2",      4,
  "Cholesterol, Fasting, High", 240,    "mg/dL",   "2",      5,
  "Cholesterol, Fasting, High", 239,    "mg/dL",   "1",      6,
  "Cholesterol, Fasting, High", 200,    "mg/dL",   "1",      7,
  "Cholesterol, Fasting, High", 199,    "mg/dL",   "0",      8,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 200,    NA,        NA,       9,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mg/dL",   NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

### < 18 years of age
### SI unit is mmol/L
### Grade 3: >= 7.77 mmol/L
### Grade 2: 5.15 -< 7.77 mmol/L
### Grade 1: 4.4 -< 5.15 mmol/L

### CV unit is mg/dL
### Grade 3: >= 300 mg/dL
### Grade 2: 200 -< 300 mg/dL
### Grade 1: 170 -< 200 mg/dL

expected_cholfilt18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 3.5,    "mmol/L",  NA,       13,
  NA_character_,                3.5,    "mmol/L",  NA,       14,
  "Cholesterol, Fasting, High", 7.77,   "mmol/L",  "3",      15,
  "Cholesterol, Fasting, High", 7.76,   "mmol/L",  "2",      16,
  "Cholesterol, Fasting, High", 5.15,   "mmol/L",  "2",      17,
  "Cholesterol, Fasting, High", 5.14,   "mmol/L",  "1",      18,
  "Cholesterol, Fasting, High", 4.4,    "mmol/L",  "1",      19,
  "Cholesterol, Fasting, High", 4.39,   "mmol/L",  "0",      20,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 3.5,    NA,        NA,       21,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mmol/L",  NA,       22,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-07")
  )

expected_cholfilt18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",                 300,    "mg/dL",  NA,       13,
  NA_character_,                300,    "mg/dL",  NA,       14,
  "Cholesterol, Fasting, High", 300,    "mg/dL",  "3",      15,
  "Cholesterol, Fasting, High", 299,    "mg/dL",  "2",      16,
  "Cholesterol, Fasting, High", 200,    "mg/dL",  "2",      17,
  "Cholesterol, Fasting, High", 199,    "mg/dL",  "1",      18,
  "Cholesterol, Fasting, High", 170,    "mg/dL",  "1",      19,
  "Cholesterol, Fasting, High", 169,    "mg/dL",  "0",      20,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 169,    NA,       NA,       21,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mg/dL",  NA,       22,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-07")
  )

## create SI records with age missing
expected_cholfinoage_daids_si <- expected_cholfige18y_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## create CV records with age missing
expected_cholfinoage_daids_cv <- expected_cholfige18y_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

## put all SI data together
expected_cholfi_daids_si <- expected_cholfinoage_daids_si %>%
  bind_rows(
    expected_cholfige18y_daids_si,
    expected_cholfilt18y_daids_si
  )

## put all CV data together
expected_cholfi_daids_cv <- expected_cholfinoage_daids_cv %>%
  bind_rows(
    expected_cholfige18y_daids_cv,
    expected_cholfilt18y_daids_cv
  )

## Test 103a: DAIDS Cholesterol, Fasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 103a: DAIDS Cholesterol, Fasting, High (SI unit)", {
  input_cholfi_daids_si <- expected_cholfi_daids_si %>%
    select(-ATOXGRH)

  actual_cholfi_daids_si <- derive_var_atoxgr_dir(
    input_cholfi_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cholfi_daids_si,
    compare = actual_cholfi_daids_si,
    keys = c("TESTNUM")
  )
})

## Test 103b: DAIDS Cholesterol, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 103b: DAIDS Cholesterol, Fasting, High (USCV unit)", {
  input_cholfi_daids_cv <- expected_cholfi_daids_cv %>%
    select(-ATOXGRH)

  actual_cholfi_daids_cv <- derive_var_atoxgr_dir(
    input_cholfi_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cholfi_daids_cv,
    compare = actual_cholfi_daids_cv,
    keys = c("TESTNUM")
  )
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

expected_ldlfige18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",         3.1,    "mmol/L",  NA,       1,
  NA_character_,        3.1,    "mmol/L",  NA,       2,
  "LDL, Fasting, High", 4.9,    "mmol/L",  "3",      3,
  "LDL, Fasting, High", 4.89,   "mmol/L",  "2",      4,
  "LDL, Fasting, High", 4.12,   "mmol/L",  "2",      5,
  "LDL, Fasting, High", 4.11,   "mmol/L",  "1",      6,
  "LDL, Fasting, High", 3.37,   "mmol/L",  "1",      7,
  "LDL, Fasting, High", 3.36,   "mmol/L",  "0",      8,
  # Unit missing cannot grade
  "LDL, Fasting, High", 3.1,    NA,        NA,       9,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mmol/L",  NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

expected_ldlfige18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,             ~AVAL,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",          190,    "mg/dL",        NA,  1,
  NA_character_,         190,    "mg/dL",        NA,  2,
  "LDL, Fasting, High",  190,    "mg/dL",       "3",  3,
  "LDL, Fasting, High",  189,    "mg/dL",       "2",  4,
  "LDL, Fasting, High",  160,    "mg/dL",       "2",  5,
  "LDL, Fasting, High",  159,    "mg/dL",       "1",  6,
  "LDL, Fasting, High",  130,    "mg/dL",       "1",  7,
  "LDL, Fasting, High",  129,    "mg/dL",       "0",  8,
  # Unit missing cannot grade
  "LDL, Fasting, High",  129,         NA,        NA,  9,
  # AVAL missing cannot grade
  "LDL, Fasting, High",  NA,     "mg/dL",        NA,  10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

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
expected_ldlfilt18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",         2.8,    "mmol/L",  NA,       11,
  NA_character_,        2.8,    "mmol/L",  NA,       12,
  "LDL, Fasting, High", 4.9,    "mmol/L",  "3",      13,
  "LDL, Fasting, High", 4.89,   "mmol/L",  "2",      14,
  "LDL, Fasting, High", 3.34,   "mmol/L",  "2",      15,
  "LDL, Fasting, High", 3.33,   "mmol/L",  "1",      16,
  "LDL, Fasting, High", 2.85,   "mmol/L",  "1",      17,
  "LDL, Fasting, High", 2.84,   "mmol/L",  "0",      18,
  # Unit missing cannot grade
  "LDL, Fasting, High", 3.5,    NA,        NA,       19,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mmol/L",  NA,       20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2020-01-07"),
    ADT = lubridate::ymd("2023-01-07")
  )
## CV data
expected_ldlfilt18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",         190,    "mg/dL",  NA,       11,
  NA_character_,        190,    "mg/dL",  NA,       12,
  "LDL, Fasting, High", 190,    "mg/dL",  "3",      13,
  "LDL, Fasting, High", 189,    "mg/dL",  "2",      14,
  "LDL, Fasting, High", 130,    "mg/dL",  "2",      15,
  "LDL, Fasting, High", 129,    "mg/dL",  "1",      16,
  "LDL, Fasting, High", 110,    "mg/dL",  "1",      17,
  "LDL, Fasting, High", 109,    "mg/dL",  "0",      18,
  # Unit missing cannot grade
  "LDL, Fasting, High", 110,    NA,       NA,       19,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mg/dL",  NA,       20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2020-01-07"),
    ADT = lubridate::ymd("2023-01-07")
  )

## create SI records with age mssing
expected_ldlfinoage_daids_si <- expected_ldlfige18y_daids_si %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    ADT = if_else(TESTNUM == 7, NA, ADT),
    BRTHDT = if_else(TESTNUM == 8, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 25, 26)
  )

## create CV records with age mssing
expected_ldlfinoage_daids_cv <- expected_ldlfige18y_daids_cv %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    ADT = if_else(TESTNUM == 7, NA, ADT),
    BRTHDT = if_else(TESTNUM == 8, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 25, 26)
  )

## create SI records with age <= 2 years
expected_ldlfile2y_daids_si <- expected_ldlfige18y_daids_si %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 7, lubridate::ymd("2021-01-07"), lubridate::ymd("2022-01-07")),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 27, 28)
  )

## create CV records with age <= 2 years
expected_ldlfile2y_daids_cv <- expected_ldlfige18y_daids_cv %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 7, lubridate::ymd("2021-01-07"), lubridate::ymd("2022-01-07")),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 27, 28)
  )

## put all SI records together
expected_ldlfi_daids_si <- expected_ldlfile2y_daids_si %>%
  bind_rows(
    expected_ldlfinoage_daids_si,
    expected_ldlfige18y_daids_si,
    expected_ldlfilt18y_daids_si
  )

## put all CV records together
expected_ldlfi_daids_cv <- expected_ldlfile2y_daids_cv %>%
  bind_rows(
    expected_ldlfinoage_daids_cv,
    expected_ldlfige18y_daids_cv,
    expected_ldlfilt18y_daids_cv
  )

## Test 104a: DAIDS LDL, Fasting, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 104a: DAIDS LDL, Fasting, High (SI unit)", {
  input_ldlfi_daids_si <- expected_ldlfi_daids_si %>%
    select(-ATOXGRH)

  actual_ldlfi_daids_si <- derive_var_atoxgr_dir(
    input_ldlfi_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ldlfi_daids_si,
    compare = actual_ldlfi_daids_si,
    keys = c("TESTNUM")
  )
})

## Test 104b: DAIDS LDL, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 104b: DAIDS LDL, Fasting, High (USCV unit)", {
  input_ldlfi_daids_cv <- expected_ldlfi_daids_cv %>%
    select(-ATOXGRH)

  actual_ldlfi_daids_cv <- derive_var_atoxgr_dir(
    input_ldlfi_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ldlfi_daids_cv,
    compare = actual_ldlfi_daids_cv,
    keys = c("TESTNUM")
  )
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
  expected_trigfi_daids_si <- tibble::tribble(
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

  input_trigfi_daids_si <- expected_trigfi_daids_si %>%
    select(-ATOXGRH)


  actual_trigfi_daids_si <- derive_var_atoxgr_dir(
    input_trigfi_daids_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigfi_daids_si,
    compare = actual_trigfi_daids_si,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 105b: DAIDS Triglycerides, Fasting, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 105b: DAIDS Triglycerides, Fasting, High (USCV unit)", {
  expected_trigfi_daids_cv <- tibble::tribble(
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

  input_trigfi_daids_cv <- expected_trigfi_daids_cv %>%
    select(-ATOXGRH)

  actual_trigfi_daids_cv <- derive_var_atoxgr_dir(
    input_trigfi_daids_cv,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigfi_daids_cv,
    compare = actual_trigfi_daids_cv,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### Magnesium, Low

### SI unit is mmol/L
### Grade 4: <0.3 mmol/L
### Grade 3: 0.3 -< 0.45 mmol/L
### Grade 2: 0.45 -< 0.6 mmol/L
### Grade 1: 0.6 -< 0.7 mmol/L

### CV unit is mg/dL
### Conversion factor 2.431 used to convert mmol/L to mg/dL

expected_magd_daids_si <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",      0.5,    "mmol/L",  NA,       1,
  NA_character_,     0.5,    "mmol/L",  NA,       2,
  "Magnesium, Low",  0.29,   "mmol/L",  "4",      3,
  "Magnesium, Low",  0.3,    "mmol/L",  "3",      4,
  "Magnesium, Low",  0.44,   "mmol/L",  "3",      5,
  "Magnesium, Low",  0.45,   "mmol/L",  "2",      6,
  "Magnesium, Low",  0.59,   "mmol/L",  "2",      7,
  "Magnesium, Low",  0.6,    "mmol/L",  "1",      8,
  "Magnesium, Low",  0.69,   "mmol/L",  "1",      9,
  "Magnesium, Low",  0.7,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Magnesium, Low",  0.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Magnesium, Low",  NA,     "mmol/L",  NA,       12,
)

## CV unit using conversions factor 2.431
expected_magd_daids_cv <- expected_magd_daids_si %>%
  mutate(
    AVAL = AVAL * 2.431,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )

## Test 106a: DAIDS Magnesium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 106a: DAIDS Magnesium, Low (SI unit)", {
  input_magd_daids_si <- expected_magd_daids_si %>%
    select(-ATOXGRL)

  actual_magd_daids_si <- derive_var_atoxgr_dir(
    input_magd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magd_daids_si,
    compare = actual_magd_daids_si,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 106b: DAIDS Magnesium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 106b: DAIDS Magnesium, Low (USCV unit)", {
  input_magd_daids_cv <- expected_magd_daids_cv %>%
    select(-ATOXGRL)

  actual_magd_daids_cv <- derive_var_atoxgr_dir(
    input_magd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magd_daids_cv,
    compare = actual_magd_daids_cv,
    keys = c("ATOXDSCL", "TESTNUM")
  )
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

expected_phosd_daids_gt14y_si <- tibble::tribble(
  ~AVAL, ~ANRLO, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  0.9,   0.8,    "MM3",    NA,       1,
  0.31,  0.8,    "mmol/L", "4",      2,
  0.32,  0.8,    "mmol/L", "3",      3,
  0.44,  0.8,    "mmol/L", "3",      4,
  0.45,  0.8,    "mmol/L", "2",      5,
  0.64,  0.8,    "mmol/L", "2",      6,
  0.65,  0.8,    "mmol/L", "1",      7,
  0.79,  0.8,    "mmol/L", "1",      8,
  0.8,   0.8,    "mmol/L", "0",      9,
  # missing ANRLO - can grade 2 - 4
  0.31,  0.8,    "mmol/L", "4",      10,
  0.32,  0.8,    "mmol/L", "3",      11,
  0.44,  0.8,    "mmol/L", "3",      12,
  0.45,  0.8,    "mmol/L", "2",      13,
  0.64,  0.8,    "mmol/L", "2",      14,
  # missing ANRLO - can grade 0 - 1
  0.65,  0.8,    "mmol/L", "1",      15,
  0.79,  0.8,    "mmol/L", "1",      16,
  0.8,   0.8,    "mmol/L", "0",      17,
  # missing AVAL
  NA,    0.8,    "mmol/L", NA,       18,
  # missing UNIT
  1,     0.8,    NA,       NA,       19,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2008-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

expected_phosd_daids_gt14y_cv <- tibble::tribble(
  ~AVAL,  ~ANRLO,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  0.9,    2.5,       "MM3",        NA,  1,
  0.9,    2.5,     "mg/dL",       "4",  2,
  1,      2.5,     "mg/dL",       "3",  3,
  1.3,    2.5,     "mg/dL",       "3",  4,
  1.4,    2.5,     "mg/dL",       "2",  5,
  1.9,    2.5,     "mg/dL",       "2",  6,
  2,      2.5,     "mg/dL",       "1",  7,
  2.49,   2.5,     "mg/dL",       "1",  8,
  2.5,    2.5,     "mg/dL",       "0",  9,
  # missing ANRLO - can grade 2 - 4
  0.9,    2.5,     "mg/dL",       "4",  10,
  1,      2.5,     "mg/dL",       "3",  11,
  1.3,    2.5,     "mg/dL",       "3",  12,
  1.4,    2.5,     "mg/dL",       "2",  13,
  1.9,    2.5,     "mg/dL",       "2",  14,
  # missing ANRLO - can grade 0 - 1
  2,      2.5,     "mg/dL",       "1",  15,
  2.49,   2.5,     "mg/dL",       "1",  16,
  2.5,    2.5,     "mg/dL",       "0",  17,
  # missing AVAL
  NA,     2.5,     "mg/dL",        NA,  18,
  # missing UNIT
  1,      2.5,          NA,        NA,  19,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2008-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

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

expected_phosd_daids_le14y_si <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.2,    "MM3",     NA,       20,
  0.47,   "mmol/L",  "4",      21,
  0.48,   "mmol/L",  "3",      22,
  0.8,    "mmol/L",  "3",      23,
  0.81,   "mmol/L",  "2",      24,
  0.96,   "mmol/L",  "2",      25,
  0.97,   "mmol/L",  "1",      26,
  1.12,   "mmol/L",  "1",      27,
  1.13,   "mmol/L",  "0",      28,
  NA,     "mmol/L",  NA,       29,
  1,      NA,        NA,       30,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2022-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

expected_phosd_daids_le14y_cv <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.49,   "MM3",     NA,       20,
  1.49,   "mg/dL",   "4",      21,
  1.5,    "mg/dL",   "3",      22,
  2.4,    "mg/dL",   "3",      23,
  2.5,    "mg/dL",   "2",      24,
  2.9,    "mg/dL",   "2",      25,
  3,      "mg/dL",   "1",      26,
  3.49,   "mg/dL",   "1",      27,
  3.5,    "mg/dL",   "0",      28,
  NA,     "mg/dL",   NA,       29,
  3.5,    NA,        NA,       30,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2022-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

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

expected_phosd_daids_lt1y_si <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.5,    "MM3",     NA,       31,
  0.47,   "mmol/L",  "4",      32,
  0.48,   "mmol/L",  "3",      33,
  0.8,    "mmol/L",  "3",      34,
  0.81,   "mmol/L",  "2",      35,
  1.12,   "mmol/L",  "2",      36,
  1.13,   "mmol/L",  "1",      37,
  1.44,   "mmol/L",  "1",      38,
  1.45,   "mmol/L",  "0",      39,
  NA,     "mmol/L",  NA,       40,
  1.5,    NA,        NA,       41,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )

expected_phosd_daids_lt1y_cv <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.49,   "MM3",     NA,       31,
  1.49,   "mg/dL",   "4",      32,
  1.5,    "mg/dL",   "3",      33,
  2.4,    "mg/dL",   "3",      34,
  2.5,    "mg/dL",   "2",      35,
  3.4,    "mg/dL",   "2",      36,
  3.5,    "mg/dL",   "1",      37,
  4.49,   "mg/dL",   "1",      38,
  4.5,    "mg/dL",   "0",      39,
  NA,     "mg/dL",   NA,       40,
  4.5,    NA,        NA,       41,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )


# Set lab date or birth date to missing for SI records
expected_phosd_daids_noage_si <- expected_phosd_daids_gt14y_si %>%
  filter(TESTNUM %in% c(8, 9)) %>%
  mutate(
    ADT = if_else(TESTNUM == 8, NA, ADT),
    BRTHDT = if_else(TESTNUM == 9, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = if_else(TESTNUM == 8, 42, 43)
  )
# Set lab date or birth date to missing for CV records
expected_phosd_daids_noage_cv <- expected_phosd_daids_gt14y_cv %>%
  filter(TESTNUM %in% c(8, 9)) %>%
  mutate(
    ADT = if_else(TESTNUM == 8, NA, ADT),
    BRTHDT = if_else(TESTNUM == 9, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = if_else(TESTNUM == 8, 42, 43)
  )

## put all SI records together
expected_phosd_daids_si <- expected_phosd_daids_gt14y_si %>%
  bind_rows(
    expected_phosd_daids_le14y_si,
    expected_phosd_daids_lt1y_si,
    expected_phosd_daids_noage_si
  )

## put all CV records together
expected_phosd_daids_cv <- expected_phosd_daids_gt14y_cv %>%
  bind_rows(
    expected_phosd_daids_le14y_cv,
    expected_phosd_daids_lt1y_cv,
    expected_phosd_daids_noage_cv
  )

## Test 107a: DAIDS Phosphate, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 107a: DAIDS Phosphate, Low (SI unit)", {
  input_phosd_daids_si <- expected_phosd_daids_si %>%
    select(-ATOXGRL)

  actual_phosd_daids_si <- derive_var_atoxgr_dir(
    input_phosd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_phosd_daids_si,
    compare = actual_phosd_daids_si,
    keys = c("TESTNUM")
  )
})

## Test 107b: DAIDS Phosphate, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 107b: DAIDS Phosphate, Low (USCV unit)", {
  input_phosd_daids_cv <- expected_phosd_daids_cv %>%
    select(-ATOXGRL)

  actual_phosd_daids_cv <- derive_var_atoxgr_dir(
    input_phosd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_phosd_daids_cv,
    compare = actual_phosd_daids_cv,
    keys = c("TESTNUM")
  )
})

### Potassium, High

### SI unit + CV is mmol/L
### Grade 4: >= 7 mmol/L
### Grade 3: 6.5 -< 7 mmol/L
### Grade 2: 6 -< 6.5 mmol/L
### Grade 1: 5.6 -< 6 mmol/L

expected_poti_daids <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",       5,      "mmol/L",  NA,       1,
  NA_character_,      5,      "mmol/L",  NA,       2,
  "Potassium, High",  7,      "mmol/L",  "4",      3,
  "Potassium, High",  6.9,    "mmol/L",  "3",      4,
  "Potassium, High",  6.5,    "mmol/L",  "3",      5,
  "Potassium, High",  6.4,    "mmol/L",  "2",      6,
  "Potassium, High",  6,      "mmol/L",  "2",      7,
  "Potassium, High",  5.9,    "mmol/L",  "1",      8,
  "Potassium, High",  5.6,    "mmol/L",  "1",      9,
  "Potassium, High",  5.5,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Potassium, High",  5,      NA,        NA,       11,
  # AVAL missing cannot grade
  "Potassium, High",  NA,     "mmol/L",  NA,       12,
)

input_poti_daids <- expected_poti_daids %>%
  select(-ATOXGRH)

## Test 108a: DAIDS Potassium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 108a: DAIDS Potassium, High (SI unit)", {
  actual_poti_daids <- derive_var_atoxgr_dir(
    input_poti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_poti_daids,
    compare = actual_poti_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 108b: DAIDS Potassium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 108b: DAIDS Potassium, High (USCV unit)", {
  actual_poti_daids <- derive_var_atoxgr_dir(
    input_poti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_poti_daids,
    compare = actual_poti_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### Potassium, Low

### SI unit + CV is mmol/L
### Grade 4: <2 mmol/L
### Grade 3: 2 -< 2.5 mmol/L
### Grade 2: 2.5 -< 3 mmol/L
### Grade 1: 3 -< 3.4 mmol/L

expected_potd_daids <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",      3,      "mmol/L",  NA,       1,
  NA_character_,     3,      "mmol/L",  NA,       2,
  "Potassium, Low",  1.9,    "mmol/L",  "4",      3,
  "Potassium, Low",  2,      "mmol/L",  "3",      4,
  "Potassium, Low",  2.4,    "mmol/L",  "3",      5,
  "Potassium, Low",  2.5,    "mmol/L",  "2",      6,
  "Potassium, Low",  2.9,    "mmol/L",  "2",      7,
  "Potassium, Low",  3,      "mmol/L",  "1",      8,
  "Potassium, Low",  3.3,    "mmol/L",  "1",      9,
  "Potassium, Low",  3.4,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Potassium, Low",  3,      NA,        NA,       11,
  # AVAL missing cannot grade
  "Potassium, Low",  NA,     "mmol/L",  NA,       12,
)

input_potd_daids <- expected_potd_daids %>%
  select(-ATOXGRL)

## Test 109a: DAIDS Potassium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 109a: DAIDS Potassium, Low (SI unit)", {
  actual_potd_daids <- derive_var_atoxgr_dir(
    input_potd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_potd_daids,
    compare = actual_potd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 109b: DAIDS Potassium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 109b: DAIDS Potassium, Low (USCV unit)", {
  actual_potd_daids <- derive_var_atoxgr_dir(
    input_potd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_potd_daids,
    compare = actual_potd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

### Sodium, High

### SI unit + CV is mmol/L
### Grade 4: >= 160 mmol/L
### Grade 3: 154 -< 160 mmol/L
### Grade 2: 150 -< 154 mmol/L
### Grade 1: 146 -< 150 mmol/L

expected_sodi_daids <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    146,    "mmol/L",  NA,       1,
  NA_character_,   146,    "mmol/L",  NA,       2,
  "Sodium, High",  160,    "mmol/L",  "4",      3,
  "Sodium, High",  159,    "mmol/L",  "3",      4,
  "Sodium, High",  154,    "mmol/L",  "3",      5,
  "Sodium, High",  153,    "mmol/L",  "2",      6,
  "Sodium, High",  150,    "mmol/L",  "2",      7,
  "Sodium, High",  149,    "mmol/L",  "1",      8,
  "Sodium, High",  146,    "mmol/L",  "1",      9,
  "Sodium, High",  145,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Sodium, High",  140,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Sodium, High",  NA,     "mmol/L",  NA,       12,
)

input_sodi_daids <- expected_sodi_daids %>%
  select(-ATOXGRH)

## Test 110a: DAIDS Sodium, High (SI unit) ----
test_that("derive_var_atoxgr_dir Test 110a: DAIDS Sodium, High (SI unit)", {
  actual_sodi_daids <- derive_var_atoxgr_dir(
    input_sodi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_sodi_daids,
    compare = actual_sodi_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 110b: DAIDS Sodium, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 110b: DAIDS Sodium, High (USCV unit)", {
  actual_sodi_daids <- derive_var_atoxgr_dir(
    input_sodi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_sodi_daids,
    compare = actual_sodi_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

### Sodium, Low

### SI unit + CV is mmol/L
### Grade 4: <= 120 mmol/L
### Grade 3: >120 -< 125 mmol/L
### Grade 2: 125 -< 130 mmol/L
### Grade 1: 130 -< 135 mmol/L

expected_sodd_daids <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",   119,    "mmol/L",  NA,       1,
  NA_character_,  119,    "mmol/L",  NA,       2,
  "Sodium, Low",  120,    "mmol/L",  "4",      3,
  "Sodium, Low",  121,    "mmol/L",  "3",      4,
  "Sodium, Low",  124,    "mmol/L",  "3",      5,
  "Sodium, Low",  125,    "mmol/L",  "2",      6,
  "Sodium, Low",  129,    "mmol/L",  "2",      7,
  "Sodium, Low",  130,    "mmol/L",  "1",      8,
  "Sodium, Low",  134,    "mmol/L",  "1",      9,
  "Sodium, Low",  135,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Sodium, Low",  140,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Sodium, Low",  NA,     "mmol/L",  NA,       12,
)

input_sodd_daids <- expected_sodd_daids %>%
  select(-ATOXGRL)

## Test 111a: DAIDS Sodium, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 111a: DAIDS Sodium, Low (SI unit)", {
  actual_sodd_daids <- derive_var_atoxgr_dir(
    input_sodd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_sodd_daids,
    compare = actual_sodd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 111b: DAIDS Sodium, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 111b: DAIDS Sodium, Low (USCV unit)", {
  actual_sodd_daids <- derive_var_atoxgr_dir(
    input_sodd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_sodd_daids,
    compare = actual_sodd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
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
  expected_urici_daids <- tibble::tribble(
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

  input_urici_daids <- expected_urici_daids %>%
    select(-ATOXGRH)

  actual_urici_daids <- derive_var_atoxgr_dir(
    input_urici_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_daids,
    compare = actual_urici_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
})

## Test 112b: DAIDS Uric Acid, High (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 112b: DAIDS Uric Acid, High (USCV unit)", {
  expected_urici_daids <- tibble::tribble(
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

  input_urici_daids <- expected_urici_daids %>%
    select(-ATOXGRH)

  actual_urici_daids <- derive_var_atoxgr_dir(
    input_urici_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_urici_daids,
    compare = actual_urici_daids,
    keys = c("ATOXDSCH", "TESTNUM")
  )
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

expected_cd4d_daids_gt5y <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.1,    "MM3",     NA,       1,
  0.09,   "10^9/L",  "4",      2,
  0.1,    "10^9/L",  "3",      3,
  0.19,   "10^9/L",  "3",      4,
  0.2,    "10^9/L",  "2",      5,
  0.29,   "10^9/L",  "2",      6,
  0.3,    "10^9/L",  "1",      7,
  0.39,   "10^9/L",  "1",      8,
  0.4,    "10^9/L",  "0",      9,
  NA,     "10^9/L",  NA,       10,
  1,      NA,        NA,       11,
) %>%
  mutate(
    ATOXDSCL = "Absolute CD4+ Count, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2029-07-01")
  )

# no criteria for age <= 5 years set grade to missing
expected_cd4d_daids_le5y <- expected_cd4d_daids_gt5y %>%
  mutate(
    ADT = lubridate::ymd("2028-07-01"),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 11
  )

# add missing ADT and BRTHDT
expected_cd4d_daids_noage <- expected_cd4d_daids_gt5y %>%
  filter(TESTNUM %in% c(5, 6)) %>%
  mutate(
    ADT = if_else(TESTNUM == 5, NA, ADT),
    BRTHDT = if_else(TESTNUM == 6, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 5, 23, 24)
  )

expected_cd4d_daids_si <- expected_cd4d_daids_gt5y %>%
  bind_rows(
    expected_cd4d_daids_le5y,
    expected_cd4d_daids_noage
  )

input_cd4d_daids_si <- expected_cd4d_daids_si %>%
  select(-ATOXGRL)

## Test 113a: DAIDS Absolute CD4 Count, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 113a: DAIDS Absolute CD4 Count, Low (SI unit)", {
  actual_cd4d_daids_si <- derive_var_atoxgr_dir(
    input_cd4d_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4d_daids_si,
    compare = actual_cd4d_daids_si,
    keys = c("TESTNUM")
  )
})


expected_cd4d_daids_cv <- expected_cd4d_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "1/uL", AVALU)
  )

input_cd4d_daids_cv <- expected_cd4d_daids_cv %>%
  select(-ATOXGRL)

## Test 113b: DAIDS Absolute CD4 Count, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 113b: DAIDS Absolute CD4 Count, Low (USCV unit)", {
  actual_cd4d_daids_cv <- derive_var_atoxgr_dir(
    input_cd4d_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4d_daids_cv,
    compare = actual_cd4d_daids_cv,
    keys = c("TESTNUM")
  )
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

expected_lymphd_daids_gt5y <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.35,   "MM3",     NA,       1,
  0.34,   "10^9/L",  "4",      2,
  0.35,   "10^9/L",  "3",      3,
  0.49,   "10^9/L",  "3",      4,
  0.5,    "10^9/L",  "2",      5,
  0.59,   "10^9/L",  "2",      6,
  0.6,    "10^9/L",  "1",      7,
  0.64,   "10^9/L",  "1",      8,
  0.65,   "10^9/L",  "0",      9,
  0.66,   "10^9/L",  "0",      10,
  NA,     "10^9/L",  NA,       11,
  1,      NA,        NA,       12,
) %>%
  mutate(
    ATOXDSCL = "Absolute Lymphocyte Count, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2029-07-01")
  )

# no criteria for age <= 5 years set grade to missing
expected_lymphd_daids_le5y <- expected_lymphd_daids_gt5y %>%
  mutate(
    ADT = lubridate::ymd("2028-07-01"),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 12
  )

# add missing ADT and BRTHDT
expected_lymphd_daids_noage <- expected_lymphd_daids_gt5y %>%
  filter(TESTNUM %in% c(5, 6)) %>%
  mutate(
    ADT = if_else(TESTNUM == 5, NA, ADT),
    BRTHDT = if_else(TESTNUM == 6, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 5, 25, 26)
  )

expected_lymphd_daids_si <- expected_lymphd_daids_gt5y %>%
  bind_rows(
    expected_lymphd_daids_le5y,
    expected_lymphd_daids_noage
  )

input_lymphd_daids_si <- expected_lymphd_daids_si %>%
  select(-ATOXGRL)


## Test 114a: DAIDS Absolute Lymphocyte Count, Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 114a: DAIDS Absolute Lymphocyte Count, Low (SI unit)", {
  actual_lymphd_daids_si <- derive_var_atoxgr_dir(
    input_lymphd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymphd_daids_si,
    compare = actual_lymphd_daids_si,
    keys = c("TESTNUM")
  )
})

expected_lymphd_daids_cv <- expected_lymphd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

input_lymphd_daids_cv <- expected_lymphd_daids_cv %>%
  select(-ATOXGRL)

## Test 114b: DAIDS Absolute Lymphocyte Count, Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 114b: DAIDS Absolute Lymphocyte Count, Low (USCV unit)", {
  actual_lymphd_daids_cv <- derive_var_atoxgr_dir(
    input_lymphd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymphd_daids_cv,
    compare = actual_lymphd_daids_cv,
    keys = c("TESTNUM")
  )
})

expected_lymphd_daids_cv2 <- expected_lymphd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

input_lymphd_daids_cv2 <- expected_lymphd_daids_cv2 %>%
  select(-ATOXGRL)

## Test 114c: DAIDS Absolute Lymph. Count, Low (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 114c: DAIDS Absolute Lymph. Count, Low (legacy USCV unit)", {
  actual_lymphd_daids_cv2 <- derive_var_atoxgr_dir(
    input_lymphd_daids_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymphd_daids_cv2,
    compare = actual_lymphd_daids_cv2,
    keys = c("TESTNUM")
  )
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

expected_ancd_daids_gt7d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.3,    "MM3",     NA,       1,
  0.399,  "10^9/L",  "4",      2,
  0.4,    "10^9/L",  "3",      3,
  0.599,  "10^9/L",  "3",      4,
  0.6,    "10^9/L",  "2",      5,
  0.799,  "10^9/L",  "2",      6,
  0.8,    "10^9/L",  "1",      7,
  1,      "10^9/L",  "1",      8,
  1.01,   "10^9/L",  "0",      9,
  NA,     "10^9/L",  NA,       10,
  1,      NA,        NA,       11,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-09")
  )

### 2 to 7 days of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 0.75 10^9/L
### Grade 3: 0.75 to < 1.0 10^9/L
### Grade 2: 1.0 to < 1.25 10^9/L
### Grade 1: 1.25 to 1.5 10^9/L

expected_ancd_daids_ge2d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.7,    "MM3",     NA,       12,
  0.749,  "10^9/L",  "4",      13,
  0.75,   "10^9/L",  "3",      14,
  0.999,  "10^9/L",  "3",      15,
  1,      "10^9/L",  "2",      16,
  1.249,  "10^9/L",  "2",      17,
  1.25,   "10^9/L",  "1",      18,
  1.5,    "10^9/L",  "1",      19,
  1.51,   "10^9/L",  "0",      20,
  NA,     "10^9/L",  NA,       21,
  1,      NA,        NA,       22,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-03")
  )

### <= 1 day of age

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: < 1.50 10^9/L
### Grade 3: 1.50 to < 3.0 10^9/L
### Grade 2: 3.0 to < 4.0 10^9/L
### Grade 1: 4.0 to 5.0 10^9/L

expected_ancd_daids_le1d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.5,    "MM3",     NA,       23,
  1.499,  "10^9/L",  "4",      24,
  1.5,    "10^9/L",  "3",      25,
  2.999,  "10^9/L",  "3",      26,
  3,      "10^9/L",  "2",      27,
  3.999,  "10^9/L",  "2",      28,
  4,      "10^9/L",  "1",      29,
  5,      "10^9/L",  "1",      30,
  5.01,   "10^9/L",  "0",      31,
  NA,     "10^9/L",  NA,       32,
  5,      NA,        NA,       33,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )

expected_ancd_daids <- expected_ancd_daids_gt7d %>%
  bind_rows(
    expected_ancd_daids_ge2d,
    expected_ancd_daids_le1d
  )

# Set lab date/birth date to missing
expected_ancd_daids_noage <- expected_ancd_daids %>%
  filter(TESTNUM %in% c(20, 31)) %>%
  mutate(
    ADT = if_else(TESTNUM == 20, NA, ADT),
    BRTHDT = if_else(TESTNUM == 31, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = case_when(
      TESTNUM == 20 ~ 34,
      TESTNUM == 31 ~ 35
    )
  )

expected_ancd_daids_si <- expected_ancd_daids %>%
  bind_rows(expected_ancd_daids_noage)


input_ancd_daids_si <- expected_ancd_daids_si %>%
  select(-ATOXGRL)


## Test 115a: DAIDS ANC Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 115a: DAIDS ANC Low (SI unit)", {
  actual_ancd_daids_si <- derive_var_atoxgr_dir(
    input_ancd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ancd_daids_si,
    compare = actual_ancd_daids_si,
    keys = c("TESTNUM")
  )
})

expected_ancd_daids_cv <- expected_ancd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )
input_ancd_daids_cv <- expected_ancd_daids_cv %>%
  select(-ATOXGRL)

## Test 115b: DAIDS ANC Low (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 115b: DAIDS ANC Low (USCV unit)", {
  actual_ancd_daids_cv <- derive_var_atoxgr_dir(
    input_ancd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ancd_daids_cv,
    compare = actual_ancd_daids_cv,
    keys = c("TESTNUM")
  )
})

expected_ancd_daids_cv2 <- expected_ancd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
input_ancd_daids_cv2 <- expected_ancd_daids_cv2 %>%
  select(-ATOXGRL)

## Test 115c: DAIDS ANC Low (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 115c: DAIDS ANC Low (legacy USCV unit)", {
  actual_ancd_daids_cv2 <- derive_var_atoxgr_dir(
    input_ancd_daids_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ancd_daids_cv2,
    compare = actual_ancd_daids_cv2,
    keys = c("TESTNUM")
  )
})


### Fibrinogen Decreased
### SI unit is g/L
### Grade 4: <0.5 g/L OR < 0.25 x LLN
### Grade 3: 0.5 to <0.75 g/L OR 0.25 to < 0.50 x LLN
### Grade 2: 0.75 to <1 g/L OR >= 0.50 to < 0.75 x LLN
### Grade 1: 1 to < 2 g/L OR 0.75 to < 1.00 x LLN

### CV unit is mg/dL (= 100 * g/L)

expected_fibd_daids_si <- tibble::tribble(
  ~ATOXDSCL,              ~AVAL, ~ANRLO, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",           2,     1,      "g/L",  NA,       1,
  NA_character_,          2,     1,      "g/L",  NA,       2,
  "Fibrinogen Decreased", 2,     1,      "g/dL", NA,       3,
  # test first half of criteria
  "Fibrinogen Decreased", 0.49,  1,      "g/L",  "4",      4,
  "Fibrinogen Decreased", 0.5,   1,      "g/L",  "3",      5,
  "Fibrinogen Decreased", 0.74,  1,      "g/L",  "3",      6,
  "Fibrinogen Decreased", 0.75,  1,      "g/L",  "2",      7,
  "Fibrinogen Decreased", 0.99,  1,      "g/L",  "2",      8,
  "Fibrinogen Decreased", 1,     1,      "g/L",  "1",      9,
  "Fibrinogen Decreased", 1.99,  1,      "g/L",  "1",      10,
  "Fibrinogen Decreased", 2,     1,      "g/L",  "0",      11,
  # test second half of criteria
  "Fibrinogen Decreased", 0.74,  3,      "g/L",  "4",      12,
  "Fibrinogen Decreased", 0.75,  3,      "g/L",  "3",      13,
  "Fibrinogen Decreased", 1.49,  3,      "g/L",  "3",      14,
  "Fibrinogen Decreased", 1.5,   3,      "g/L",  "2",      15,
  "Fibrinogen Decreased", 2.24,  3,      "g/L",  "2",      16,
  "Fibrinogen Decreased", 2.25,  3,      "g/L",  "1",      17,
  "Fibrinogen Decreased", 2.99,  3,      "g/L",  "1",      18,
  "Fibrinogen Decreased", 3,     3,      "g/L",  "0",      19,
  # TEST for missing values
  "Fibrinogen Decreased", 0.49,  NA,     "g/L",  "4",      20,
  "Fibrinogen Decreased", 0.5,   NA,     "g/L",  NA,       21,
  "Fibrinogen Decreased", 2,     1,      NA,     NA,       22,
  "Fibrinogen Decreased", NA,    1,      "g/L",  NA,       23,
)

input_fibd_daids_si <- expected_fibd_daids_si %>%
  select(-ATOXGRL)

## Test 116a: DAIDS Fibrinogen Decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 116a: DAIDS Fibrinogen Decreased (SI unit)", {
  actual_fibd_daids_si <- derive_var_atoxgr_dir(
    input_fibd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_fibd_daids_si,
    compare = actual_fibd_daids_si,
    keys = c("TESTNUM")
  )
})

## CV unit is mg/dL which is 100 * SI unit of g/L
expected_fibd_daids_cv <- expected_fibd_daids_si %>%
  mutate(
    AVAL = 100 * AVAL,
    ANRLO = 100 * ANRLO,
    AVALU = if_else(str_to_upper(AVALU) == "G/L", "mg/dL", NA_character_)
  )

input_fibd_daids_cv <- expected_fibd_daids_cv %>%
  select(-ATOXGRL)

## Test 116b: DAIDS Fibrinogen Decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 116b: DAIDS Fibrinogen Decreased (USCV unit)", {
  actual_fibd_daids_cv <- derive_var_atoxgr_dir(
    input_fibd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_fibd_daids_cv,
    compare = actual_fibd_daids_cv,
    keys = c("TESTNUM")
  )
})

### Hemoglobin, Low
### SI unit is g/L
### CV unit g/dL ( = g/L / 10)

### >= 13 years of age (male only)

### Grade 4: < 70 g/L
### Grade 3: 70 to < 90 g/L
### Grade 2: 90 to < 100 g/L
### Grade 1: 100 to 109 g/L

expected_hgbd_daids_ge13ym <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~SEX, ~TESTNUM,
  "Hemoglobin, Low", 69,    "MM3",  NA,       "M",  1,
  "Hemoglobin, Low", 69,    "g/L",  "4",      "M",  2,
  "Hemoglobin, Low", 70,    "g/L",  "3",      "M",  3,
  "Hemoglobin, Low", 89,    "g/L",  "3",      "M",  4,
  "Hemoglobin, Low", 90,    "g/L",  "2",      "M",  5,
  "Hemoglobin, Low", 99,    "g/L",  "2",      "M",  6,
  "Hemoglobin, Low", 100,   "g/L",  "1",      "M",  7,
  "Hemoglobin, Low", 109,   "g/L",  "1",      "M",  8,
  "Hemoglobin, Low", 110,   "g/L",  "0",      "M",  9,
  "Hemoglobin, Low", NA,    "g/L",  NA,       "M",  10,
  "Hemoglobin, Low", 110,   NA,     NA,       "M",  11,
  "Hemoglobin, Low", 110,   "g/L",  NA,       NA,   12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

### >= 13 years of age (female only)

### Grade 4: < 65 g/L
### Grade 3: 65 to < 85 g/L
### Grade 2: 85 to < 95 g/L
### Grade 1: 95 to 104 g/L

expected_hgbd_daids_ge13yf <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~SEX, ~TESTNUM,
  "Hemoglobin, Low", 64,    "MM3",  NA,       "F",  13,
  "Hemoglobin, Low", 64,    "g/L",  "4",      "F",  14,
  "Hemoglobin, Low", 65,    "g/L",  "3",      "F",  15,
  "Hemoglobin, Low", 84,    "g/L",  "3",      "F",  16,
  "Hemoglobin, Low", 85,    "g/L",  "2",      "F",  17,
  "Hemoglobin, Low", 94,    "g/L",  "2",      "F",  18,
  "Hemoglobin, Low", 95,    "g/L",  "1",      "F",  19,
  "Hemoglobin, Low", 104,   "g/L",  "1",      "F",  20,
  "Hemoglobin, Low", 105,   "g/L",  "0",      "F",  21,
  "Hemoglobin, Low", NA,    "g/L",  NA,       "F",  22,
  "Hemoglobin, Low", 110,   NA,     NA,       "F",  23,
  "Hemoglobin, Low", 110,   "g/L",  NA,       NA,   24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )


### 57 days to < 13 years of age (male and female)

### Grade 4: < 65 g/L
### Grade 3: 65 to < 85 g/L
### Grade 2: 85 to < 95 g/L
### Grade 1: 95 to 104 g/L

expected_hgbd_daids_lt13y <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 64,    "MM3",  NA,       25,
  "Hemoglobin, Low", 64,    "g/L",  "4",      26,
  "Hemoglobin, Low", 65,    "g/L",  "3",      27,
  "Hemoglobin, Low", 84,    "g/L",  "3",      28,
  "Hemoglobin, Low", 85,    "g/L",  "2",      29,
  "Hemoglobin, Low", 94,    "g/L",  "2",      30,
  "Hemoglobin, Low", 95,    "g/L",  "1",      31,
  "Hemoglobin, Low", 104,   "g/L",  "1",      32,
  "Hemoglobin, Low", 105,   "g/L",  "0",      33,
  "Hemoglobin, Low", NA,    "g/L",  NA,       34,
  "Hemoglobin, Low", 110,   NA,     NA,       35,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-06-30")
  )

### 36 to <= 56 days of age (male and female)

### Grade 4: < 60 g/L
### Grade 3: 60 to < 70 g/L
### Grade 2: 70 to < 85 g/L
### Grade 1: 85 to 96 g/L

expected_hgbd_daids_le56d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 59,    "MM3",  NA,       36,
  "Hemoglobin, Low", 59,    "g/L",  "4",      37,
  "Hemoglobin, Low", 60,    "g/L",  "3",      38,
  "Hemoglobin, Low", 69,    "g/L",  "3",      39,
  "Hemoglobin, Low", 70,    "g/L",  "2",      40,
  "Hemoglobin, Low", 84,    "g/L",  "2",      41,
  "Hemoglobin, Low", 85,    "g/L",  "1",      42,
  "Hemoglobin, Low", 96,    "g/L",  "1",      43,
  "Hemoglobin, Low", 97,    "g/L",  "0",      44,
  "Hemoglobin, Low", NA,    "g/L",  NA,       45,
  "Hemoglobin, Low", 110,   NA,     NA,       46,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-08-26")
  )


### 22 to <= 35 days of age (male and female)

### Grade 4: < 67 g/L
### Grade 3: 67 to < 80 g/L
### Grade 2: 80 to < 95 g/L
### Grade 1: 95 to 110 g/L

expected_hgbd_daids_le35d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 66,    "MM3",  NA,       47,
  "Hemoglobin, Low", 66,    "g/L",  "4",      48,
  "Hemoglobin, Low", 67,    "g/L",  "3",      49,
  "Hemoglobin, Low", 79,    "g/L",  "3",      50,
  "Hemoglobin, Low", 80,    "g/L",  "2",      51,
  "Hemoglobin, Low", 94,    "g/L",  "2",      52,
  "Hemoglobin, Low", 95,    "g/L",  "1",      53,
  "Hemoglobin, Low", 110,   "g/L",  "1",      54,
  "Hemoglobin, Low", 111,   "g/L",  "0",      55,
  "Hemoglobin, Low", NA,    "g/L",  NA,       56,
  "Hemoglobin, Low", 110,   NA,     NA,       57,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-08-05")
  )


### 8 to <= 21 days of age (male and female)

### Grade 4: < 80 g/L
### Grade 3: 80 to < 90 g/L
### Grade 2: 90 to < 110 g/L
### Grade 1: 110 to 130 g/L

expected_hgbd_daids_le21d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 79,    "MM3",  NA,       58,
  "Hemoglobin, Low", 79,    "g/L",  "4",      59,
  "Hemoglobin, Low", 80,    "g/L",  "3",      60,
  "Hemoglobin, Low", 89,    "g/L",  "3",      61,
  "Hemoglobin, Low", 90,    "g/L",  "2",      62,
  "Hemoglobin, Low", 109,   "g/L",  "2",      63,
  "Hemoglobin, Low", 110,   "g/L",  "1",      64,
  "Hemoglobin, Low", 130,   "g/L",  "1",      65,
  "Hemoglobin, Low", 131,   "g/L",  "0",      66,
  "Hemoglobin, Low", NA,    "g/L",  NA,       67,
  "Hemoglobin, Low", 110,   NA,     NA,       68,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-22")
  )

### <= 7 days of age (male and female)

### Grade 4: < 90 g/L
### Grade 3: 90 to < 100 g/L
### Grade 2: 100 to < 130 g/L
### Grade 1: 130 to 140 g/L

expected_hgbd_daids_le7d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 89,    "MM3",  NA,       69,
  "Hemoglobin, Low", 89,    "g/L",  "4",      70,
  "Hemoglobin, Low", 90,    "g/L",  "3",      71,
  "Hemoglobin, Low", 99,    "g/L",  "3",      72,
  "Hemoglobin, Low", 100,   "g/L",  "2",      73,
  "Hemoglobin, Low", 129,   "g/L",  "2",      74,
  "Hemoglobin, Low", 130,   "g/L",  "1",      75,
  "Hemoglobin, Low", 140,   "g/L",  "1",      76,
  "Hemoglobin, Low", 141,   "g/L",  "0",      77,
  "Hemoglobin, Low", NA,    "g/L",  NA,       78,
  "Hemoglobin, Low", 110,   NA,     NA,       79,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-08")
  )

expected_hgbd_daids <- expected_hgbd_daids_ge13ym %>%
  bind_rows(
    expected_hgbd_daids_ge13yf,
    expected_hgbd_daids_lt13y,
    expected_hgbd_daids_le56d,
    expected_hgbd_daids_le35d,
    expected_hgbd_daids_le21d,
    expected_hgbd_daids_le7d
  )

# Set lab date to missing fo each type, ie SEX is M, F or missing
expected_hgbd_daids_noage <- expected_hgbd_daids %>%
  filter(TESTNUM %in% c(5, 17, 29)) %>%
  mutate(
    ADT = NA,
    ATOXGRL = NA,
    TESTNUM = case_when(
      TESTNUM == 5 ~ 80,
      TESTNUM == 17 ~ 81,
      TESTNUM == 29 ~ 82
    )
  )

expected_hgbd_daids <- expected_hgbd_daids %>%
  bind_rows(expected_hgbd_daids_noage)


input_hgbd_daids <- expected_hgbd_daids %>%
  select(-ATOXGRL)

## Test 117a: DAIDS HGB Low (SI unit) ----
test_that("derive_var_atoxgr_dir Test 117: DAIDS HGB Low (SI unit)", {
  actual_hgbd_daids <- derive_var_atoxgr_dir(
    input_hgbd_daids,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbd_daids,
    compare = actual_hgbd_daids,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})

## Test 117b: DAIDS HGB Low (CV unit) ----
test_that("derive_var_atoxgr_dir Test 117b: DAIDS HGB Low (CV unit)", {
  expected_hgbd_daids_cv <- expected_hgbd_daids %>%
    mutate(
      AVAL = AVAL / 10,
      AVALU = if_else(!is.na(AVALU) & str_to_upper(AVALU) == "G/L", "g/dL", NA_character_)
    )

  input_hgbd_daids_cv <- expected_hgbd_daids_cv %>%
    select(-ATOXGRL)


  actual_hgbd_daids_cv <- derive_var_atoxgr_dir(
    input_hgbd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbd_daids_cv,
    compare = actual_hgbd_daids_cv,
    keys = c("ATOXDSCL", "TESTNUM")
  )
})


### INR, high

### Grade 4: >=3 x ULN
### Grade 3: 2 to <3 x ULN
### Grade 2: 1.5 to < 2 x ULN
### Grade 1: 1.1 to < 1.5 x ULN

## Test 118: DAIDS INR High ----
test_that("derive_var_atoxgr_dir Test 118: DAIDS INR High", {
  expected_inri_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",  80,     80,     NA_character_,  NA,
    NA_character_, 60,     80,     NA_character_,  NA,
    "INR, High",   240,    80,     NA_character_,  "4",
    "INR, High",   239,    80,     NA_character_,  "3",
    "INR, High",   160,    80,     NA_character_,  "3",
    "INR, High",   159,    80,     NA_character_,  "2",
    "INR, High",   120,    80,     NA_character_,  "2",
    "INR, High",   119,    80,     NA_character_,  "1",
    "INR, High",   88,     80,     NA_character_,  "1",
    "INR, High",   87,     80,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "INR, High",   100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "INR, High",   NA,     80,     NA_character_,  NA,
  )

  input_inri_daids <- expected_inri_daids %>%
    select(-ATOXGRH)

  actual_inri_daids <- derive_var_atoxgr_dir(
    input_inri_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_inri_daids,
    compare = actual_inri_daids,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})


### Methemoglobin
### SI + CV unit is %

### Grade 4: >=20.0%
### Grade 3: 15 to < 20%
### Grade 2: 10 to < 15%
### Grade 1: 5 to <10%

expected_methi_daids <- tibble::tribble(
  ~ATOXDSCH,         ~AVAL, ~AVALU,  ~ATOXGRH,
  "Not a term",      20,    "%",     NA,
  NA_character_,     20,    "%",     NA,
  "Methemoglobin",   20,    "%",     "4",
  "Methemoglobin",   19,    "%",     "3",
  "Methemoglobin",   15,    "%",     "3",
  "Methemoglobin",   14,    "%",     "2",
  "Methemoglobin",   10,    "%",     "2",
  "Methemoglobin",   9,     "%",     "1",
  "Methemoglobin",   5,     "%",     "1",
  "Methemoglobin",   4,     "%",     "0",
  # Unit wrong - cannot grade
  "Methemoglobin",   100,   NA,      NA,
  # AVAL missing cannot grade
  "Methemoglobin",   NA,    "%",     NA,
)

input_methi_daids <- expected_methi_daids %>%
  select(-ATOXGRH)

## Test 119a: DAIDS Methemoglobin (SI unit) ----
test_that("derive_var_atoxgr_dir Test 119a: DAIDS Methemoglobin (SI unit)", {
  actual_methi_daids <- derive_var_atoxgr_dir(
    input_methi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_methi_daids,
    compare = actual_methi_daids,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

## Test 119b: DAIDS Methemoglobin (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 119b: DAIDS Methemoglobin (USCV unit)", {
  actual_methi_daids <- derive_var_atoxgr_dir(
    input_methi_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_methi_daids,
    compare = actual_methi_daids,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

### PTT, high

### Grade 4: >=3 x ULN
### Grade 3: 2.33 to <3 x ULN
### Grade 2: 1.66 to < 2.33 x ULN
### Grade 1: 1.1 to < 1.66 x ULN

## Test 120: DAIDS PTT High ----
test_that("derive_var_atoxgr_dir Test 120: DAIDS PTT High", {
  expected_ptti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL, ~ANRHI, ~AVALU,        ~ATOXGRH,
    "Not a term",  80,    80,     NA_character_, NA,
    NA_character_, 60,    80,     NA_character_, NA,
    "PTT, High",   240,   80,     NA_character_, "4",
    "PTT, High",   239,   80,     NA_character_, "3",
    "PTT, High",   186.4, 80,     NA_character_, "3",
    "PTT, High",   186.3, 80,     NA_character_, "2",
    "PTT, High",   132.8, 80,     NA_character_, "2",
    "PTT, High",   132.7, 80,     NA_character_, "1",
    "PTT, High",   88,    80,     NA_character_, "1",
    "PTT, High",   87,    80,     NA_character_, "0",
    # ANRHI missing - cannot grade
    "PTT, High",   100,   NA,     NA_character_, NA,
    # AVAL missing cannot grade
    "PTT, High",   NA,    80,     NA_character_, NA,
  )

  input_ptti_daids <- expected_ptti_daids %>%
    select(-ATOXGRH)

  actual_ptti_daids <- derive_var_atoxgr_dir(
    input_ptti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ptti_daids,
    compare = actual_ptti_daids,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})


### Platelets, Decreased
### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <25 x 10e9 /L
### Grade 3: 25 to <50  x 10e9 /L
### Grade 2: 50 to <100 - x 10e9
### Grade 1: 100 - 125 x 10e9 /L

expected_plated_daids_si <- tibble::tribble(
  ~ATOXDSCL,              ~AVAL, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",           126,   "10^9/L",  NA,       1,
  NA_character_,          120,   "10^9/L",  NA,       2,
  "Platelets, Decreased", 115,   "MM3",     NA,       3,
  "Platelets, Decreased", 24.9,  "10^9/L",  "4",      4,
  "Platelets, Decreased", 25,    "10^9/L",  "3",      5,
  "Platelets, Decreased", 49.9,  "10^9/L",  "3",      6,
  "Platelets, Decreased", 50,    "10^9/L",  "2",      7,
  "Platelets, Decreased", 99.9,  "10^9/L",  "2",      8,
  "Platelets, Decreased", 100,   "10^9/L",  "1",      9,
  "Platelets, Decreased", 124.9, "10^9/L",  "1",      10,
  "Platelets, Decreased", 125,   "10^9/L",  "0",      11,
)

input_plated_daids_si <- expected_plated_daids_si %>%
  select(-ATOXGRL)

## Test 121a: DAIDS Platelets decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 121a: DAIDS Platelets decreased (SI unit)", {
  actual_plated_daids_si <- derive_var_atoxgr_dir(
    input_plated_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plated_daids_si,
    compare = actual_plated_daids_si,
    keys = c("TESTNUM")
  )
})

expected_plated_daids_cv <- expected_plated_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

input_plated_daids_cv <- expected_plated_daids_cv %>%
  select(-ATOXGRL)

## Test 121b: DAIDS Platelets decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 121b: DAIDS Platelets decreased (USCV unit)", {
  actual_plated_daids_cv <- derive_var_atoxgr_dir(
    input_plated_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plated_daids_cv,
    compare = actual_plated_daids_cv,
    keys = c("TESTNUM")
  )
})

expected_plated_daids_cv2 <- expected_plated_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

input_plated_daids_cv2 <- expected_plated_daids_cv2 %>%
  select(-ATOXGRL)

## Test 121c: DAIDS Platelets decreased (legacy USCV unit) ----
test_that("derive_var_atoxgr_dir Test 121c: DAIDS Platelets decreased (legacy USCV unit)", {
  actual_plated_daids_cv2 <- derive_var_atoxgr_dir(
    input_plated_daids_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plated_daids_cv2,
    compare = actual_plated_daids_cv2,
    keys = c("TESTNUM")
  )
})

### PT, high

### Grade 4: >=3 x ULN
### Grade 3: 1.5 - <3 x ULN
### Grade 2: 1.25 - <1.5 x ULN
### Grade 1: 1.1 - <1.25 x ULN

## Test 122: DAIDS PT High ----
test_that("derive_var_atoxgr_dir Test 122: DAIDS PT High", {
  expected_pti_daids <- tibble::tribble(
    ~ATOXDSCH,     ~AVAL,  ~ANRHI,  ~AVALU,         ~ATOXGRH,
    "Not a term",  80,     100,     NA_character_,  NA,
    NA_character_, 60,     100,     NA_character_,  NA,
    "PT, High",    300,    100,     NA_character_,  "4",
    "PT, High",    299,    100,     NA_character_,  "3",
    "PT, High",    150,    100,     NA_character_,  "3",
    "PT, High",    149,    100,     NA_character_,  "2",
    "PT, High",    125,    100,     NA_character_,  "2",
    "PT, High",    124,    100,     NA_character_,  "1",
    "PT, High",    110,    100,     NA_character_,  "1",
    "PT, High",    109,    100,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "PT, High",    100,    NA,      NA_character_,  NA,
    # AVAL missing cannot grade
    "PT, High",    NA,     100,     NA_character_,  NA,
  )

  input_pti_daids <- expected_pti_daids %>%
    select(-ATOXGRH)

  actual_pti_daids <- derive_var_atoxgr_dir(
    input_pti_daids,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_pti_daids,
    compare = actual_pti_daids,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### White blood cell decreased (> 7 days of age)

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <1 x 10e9/L
### Grade 3: 1 to 1.499 x 10e9/L
### Grade 2: 1.5 to 1.999 x 10e9/L
### Grade 1: 2 to 2.499 x 10e9/L

expected_wbcd_daids_gt7d <- tibble::tribble(
  ~ATOXDSCL,        ~AVAL, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",     1,     "10^9/L",  NA,       1,
  NA_character_,    2,     "10^9/L",  NA,       2,
  "WBC, Decreased", 0.9,   "MM3",     NA,       3,
  "WBC, Decreased", 0.9,   "10^9/L",  "4",      4,
  "WBC, Decreased", 1,     "10^9/L",  "3",      5,
  "WBC, Decreased", 1.49,  "10^9/L",  "3",      6,
  "WBC, Decreased", 1.5,   "10^9/L",  "2",      7,
  "WBC, Decreased", 1.99,  "10^9/L",  "2",      8,
  "WBC, Decreased", 2,     "10^9/L",  "1",      9,
  "WBC, Decreased", 2.49,  "10^9/L",  "1",      10,
  "WBC, Decreased", 2.5,   "10^9/L",  "0",      11,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-09")
  )

### White blood cell decreased (<= 7 days of age)

### Si unit is 10^9/L
### CV unit is 10^3/uL ( = 10^9/L)
### Legacy CV unit is 10^3/mL (1000 * 10^9/L)
### Grade 4: <2.500 x 10e9/L
### Grade 3: 2.5 to 3.999 x 10e9/L
### Grade 2:   4 to 5.499 x 10e9/L
### Grade 1: 5.5 to 6.999 x 10e9/L

expected_wbcd_daids_le7d <- tibble::tribble(
  ~ATOXDSCL,        ~AVAL, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "WBC, Decreased", 2.49,  "MM3",    NA,       12,
  "WBC, Decreased", 2.49,  "10^9/L", "4",      13,
  "WBC, Decreased", 2.5,   "10^9/L", "3",      14,
  "WBC, Decreased", 3.99,  "10^9/L", "3",      15,
  "WBC, Decreased", 4,     "10^9/L", "2",      16,
  "WBC, Decreased", 5.49,  "10^9/L", "2",      17,
  "WBC, Decreased", 5.5,   "10^9/L", "1",      18,
  "WBC, Decreased", 6.99,  "10^9/L", "1",      19,
  "WBC, Decreased", 7,     "10^9/L", "0",      20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-08")
  )

expected_wbcd_daids_noage <- expected_wbcd_daids_gt7d %>%
  filter(TESTNUM %in% c(10, 11)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ADT = if_else(TESTNUM == 11, NA, ADT),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 11
  )

expected_wbcd_daids_si <- expected_wbcd_daids_gt7d %>%
  bind_rows(
    expected_wbcd_daids_le7d,
    expected_wbcd_daids_noage
  )

input_wbcd_daids_si <- expected_wbcd_daids_si %>%
  select(-ATOXGRL)

## Test 123a: DAIDS White blood cell decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 123a: DAIDS White blood cell decreased (SI unit)", {
  actual_wbcd_daids_si <- derive_var_atoxgr_dir(
    input_wbcd_daids_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_daids_si,
    compare = actual_wbcd_daids_si,
    keys = c("TESTNUM")
  )
})

expected_wbcd_daids_cv <- expected_wbcd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

input_wbcd_daids_cv <- expected_wbcd_daids_cv %>%
  select(-ATOXGRL)

## Test 123b: DAIDS White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 123b: DAIDS White blood cell decreased (USCV unit)", {
  actual_wbcd_daids_cv <- derive_var_atoxgr_dir(
    input_wbcd_daids_cv,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_daids_cv,
    compare = actual_wbcd_daids_cv,
    keys = c("TESTNUM")
  )
})

expected_wbcd_daids_cv2 <- expected_wbcd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

input_wbcd_daids_cv2 <- expected_wbcd_daids_cv2 %>%
  select(-ATOXGRL)

## Test 123c: DAIDS White blood cell decreased (USCV unit) ----
test_that("derive_var_atoxgr_dir Test 123c: DAIDS White blood cell decreased (USCV unit)", {
  actual_wbcd_daids_cv2 <- derive_var_atoxgr_dir(
    input_wbcd_daids_cv2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_daids_uscv,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_wbcd_daids_cv2,
    compare = actual_wbcd_daids_cv2,
    keys = c("TESTNUM")
  )
})

# NCICTCAEv6 new grading----

### Blood lactate dehydrogenase increased

### Grade 1: >ULN

## Test 124: CTCAEv6 LDH (SI unit) ----
test_that("derive_var_atoxgr_dir Test 124: CTCAEv6 LDH increased (SI unit)", {
  expected_ldhi_si <- tibble::tribble(
    ~ATOXDSCH,                               ~AVAL, ~ANRHI,      ~ATOXGRH,  ~TESTNUM,
    "Not a term",                              1.1,      1, NA_character_,         1,
    NA_character_,                             1.1,      1, NA_character_,         2,
    "Blood lactate dehydrogenase increased",   1.1,      1,           "1",         3,
    "Blood lactate dehydrogenase increased",     1,      1,           "0",         4,
    "Blood lactate dehydrogenase increased",    NA,      1, NA_character_,         5,
    "Blood lactate dehydrogenase increased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  input_ldhi_si <- expected_ldhi_si %>%
    select(-ATOXGRH)

  actual_ldhi_si <- derive_var_atoxgr_dir(
    input_ldhi_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ldhi_si,
    compare = actual_ldhi_si,
    keys = c("TESTNUM")
  )
})

## Test 125: CTCAEv6 HDL decreased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 125: CTCAEv6 HDL decreased (SI unit)", {
  expected_hdld_si <- tibble::tribble(
    ~ATOXDSCL,       ~AVAL, ~ANRLO,      ~ATOXGRL,  ~TESTNUM,
    "Not a term",      0.9,      1, NA_character_,         1,
    NA_character_,     0.9,      1, NA_character_,         2,
    "HDL decreased",   0.9,      1,           "1",         3,
    "HDL decreased",     1,      1,           "0",         4,
    "HDL decreased",    NA,      1, NA_character_,         5,
    "HDL decreased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  input_hdld_si <- expected_hdld_si %>%
    select(-ATOXGRL)

  actual_hdld_si <- derive_var_atoxgr_dir(
    input_hdld_si,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hdld_si,
    compare = actual_hdld_si,
    keys = c("TESTNUM")
  )
})

## Test 126: CTCAEv6 LDL increased (SI unit) ----
test_that("derive_var_atoxgr_dir Test 126: CTCAEv6 LDL increased (SI unit)", {
  expected_ldli_si <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRHI,      ~ATOXGRH,  ~TESTNUM,
    "Not a term",      1.1,      1, NA_character_,         1,
    NA_character_,     1.1,      1, NA_character_,         2,
    "LDL increased",   1.1,      1,           "1",         3,
    "LDL increased",     1,      1,           "0",         4,
    "LDL increased",    NA,      1, NA_character_,         5,
    "LDL increased",     1,     NA, NA_character_,         6,
  ) %>%
    mutate(AVALU = NA_character_)

  input_ldli_si <- expected_ldli_si %>%
    select(-ATOXGRH)

  actual_ldli_si <- derive_var_atoxgr_dir(
    input_ldli_si,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ldli_si,
    compare = actual_ldli_si,
    keys = c("TESTNUM")
  )
})

### Acidosis
### GRADE_3: pH is <7.3
### GRADE_1: pH is <normal, but >=7.3

## Test 127: CTCAEv6 Acidosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 127: CTCAEv6 Acidosis (SI unit)", {
  expected_acido <- tibble::tribble(
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

  input_acido <- expected_acido %>%
    select(-ATOXGRL)

  actual_acido <- derive_var_atoxgr_dir(
    input_acido,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_acido,
    compare = actual_acido,
    keys = c("TESTNUM")
  )
})


### Acidosis
### GRADE_3: pH is >7.5
### GRADE_1: pH is >normal, but <=7.5

## Test 128: CTCAEv6 Alkalosis (SI unit) ----
test_that("derive_var_atoxgr_dir Test 128: CTCAEv6 Alkalosis (SI unit)", {
  expected_alkalo <- tibble::tribble(
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

  input_alkalo <- expected_alkalo %>%
    select(-ATOXGRH)

  actual_alkalo <- derive_var_atoxgr_dir(
    input_alkalo,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkalo,
    compare = actual_alkalo,
    keys = c("TESTNUM")
  )
})

### Creatinine clearance decreased
### GRADE_3: 10 - 25 ml/min
### GRADE_2: 26 - 49 ml/min

## Test 129: CTCAEv6 Creat. clear. dec. (SI unit) ----
test_that("derive_var_atoxgr_dir Test 129: CTCAEv6 Creat. clear. dec. (SI unit)", {
  expected_crcl_d <- tibble::tribble(
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

  input_crcl_d <- expected_crcl_d %>%
    select(-ATOXGRL)

  actual_crcl_d <- derive_var_atoxgr_dir(
    input_crcl_d,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv6,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_crcl_d,
    compare = actual_crcl_d,
    keys = c("TESTNUM")
  )
})

## Test when deprecated abnormal_indicator used - should map to high_indicator
## Test 130: CTCAEv6  Blood bilirubin increased ----
test_that("derive_var_atoxgr_dir Test 130: CTCAEv6  Blood bilirubin increased", {
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
    base = expected_bili_ctcv6,
    compare = actual_bili_ctcv6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU")
  )
})

## Test when high_indicator not defined - should map to high_indicator
## Test 131: CTCAEv6  Blood bili incr. high_indicator not defined ----
test_that("derive_var_atoxgr_dir Test 131: CTCAEv6  Blood bili incr. high_indicator not defined", {
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
  expected_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~BNRIND,        ~ATOXGRH,
    "Creatinine increased",  241,    40,    40,     NA_character_,  "NORMAL",       "4",
  )

  input_creatn <- expected_creatn %>%
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
  atoxgr_criteria_ctcv6 <- bind_rows(
    admiral::atoxgr_criteria_ctcv6,
    admiral::atoxgr_criteria_ctcv6 %>%
      filter(TERM == "Creatinine increased") %>%
      mutate(GRADE_CRITERIA_CODE = paste0(GRADE_CRITERIA_CODE, ","))
  )

  expected_creatn <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~BNRIND,        ~ATOXGRH,
    "Creatinine increased",  241,    40,    40,     NA_character_,  "NORMAL",       "4",
  )

  input_creatn <- expected_creatn %>%
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
