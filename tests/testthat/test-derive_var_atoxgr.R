
# ---- derive_var_atoxgr, test 1: ATOXGR cannot be graded ----
test_that("derive_var_atoxgr, test 1: ATOXGR cannot be graded", {
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

# ---- derive_var_atoxgr, test 2: derive_var_atoxgr, ATOXGR = 0 (normal) ----
test_that("derive_var_atoxgr, test 2: derive_var_atoxgr, ATOXGR = 0 (normal)", {
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

# ---- derive_var_atoxgr, test 3: ATOXGR > 0 (HYPER) ----
test_that("derive_var_atoxgr, test 3: ATOXGR > 0 (HYPER)", {
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

# ---- derive_var_atoxgr, test 4: ATOXGR < 0 (HYPO) ----
test_that("derive_var_atoxgr, test 4: ATOXGR < 0 (HYPO)", {
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


expected_anemia <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL,
  "Not a term",   80,     120,    200,    "G/L",  NA,
  NA_character_,  60,     50,     100,    "G/L",  NA,
  "ANEMIA",       79,     140,    NA,     "G/L",  "3",
  "ANEMIA",       80,     140,    NA,     "G/L",  "2",
  "Anemia",       99,     140,    NA,     "G/L",  "2",
  "Anemia",       100,    140,    NA,     "G/L",  "1",
  # wrong UNIT - GRADE should be missing
  "anemia",       100,    140,    NA,     "G/dL", NA,
  "Anemia",       139,    140,    NA,     "G/L",  "1",
  "ANEMIA",       140,    140,    NA,     "g/L",  "0",
  # ANRLO missing AVAL not grade 2 or 3 - cannot grade
  "Anemia",       140,    NA,     NA,     "G/L",  NA,
  "Anemia",       139,    NA,     NA,     "G/L",  NA,
  "Anemia",       100,    NA,     NA,     "G/L",  NA,
  # ANRLO missing but AVAL satisfies grade 2
  "Anemia",       99,     NA,     NA,     "G/L",  "2",
  # ANRLO missing but AVAL satisfies grade 3
  "Anemia",       79,     NA,     NA,     "G/L",  "3",
  # Unit missing cannot grade
  "Anemia",       140,    140,    NA,     NA,     NA,
  # AVAL missing cannot grade
  "Anemia",       NA,     140,    NA,     "G/L",  NA,
)

input_anemia <- expected_anemia %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 5: CTCAEv4 Anemia ----
test_that("derive_var_atoxgr, test 5: CTCAEv4 Anemia", {
  actual_anemia_ctcv4 <- derive_var_atoxgr_dir(
    input_anemia,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia,
    compare = actual_anemia_ctcv4,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 6: CTCAEv5 Anemia ----
test_that("derive_var_atoxgr, test 6: CTCAEv5 Anemia", {
  actual_anemia_ctcv5 <- derive_var_atoxgr_dir(
    input_anemia,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_anemia,
    compare = actual_anemia_ctcv5,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})


### Leukocytosis
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: >100,000/mm3

expected_leukocytosis <- tibble::tribble(
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
input_leukocytosis <- expected_leukocytosis %>%
  select(-ATOXGRH)


# ---- derive_var_atoxgr, test 7: CTCAEv4 Leukocytosis ----
test_that("derive_var_atoxgr, test 7: CTCAEv4 Leukocytosis", {
  actual_leukocytosis <- derive_var_atoxgr_dir(
    input_leukocytosis,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis,
    compare = actual_leukocytosis,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 8: CTCAEv5 Leukocytosis ----
test_that("derive_var_atoxgr, test 8: CTCAEv5 Leukocytosis", {
  actual_leukocytosis <- derive_var_atoxgr_dir(
    input_leukocytosis,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_leukocytosis,
    compare = actual_leukocytosis,
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

# ---- derive_var_atoxgr, test 9: CTCAEv4 Activated partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr, test 9: CTCAEv4 Activated partial thromboplastin time prolonged", {
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

# ---- derive_var_atoxgr, test 10: CTCAEv5 Activated partial thromboplastin time prolonged ----
test_that("derive_var_atoxgr, test 10: CTCAEv5 Activated partial thromboplastin time prolonged", {
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

# ---- derive_var_atoxgr, test 11: CTCAEv4 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr, test 11: CTCAEv4 Alanine aminotransferase increased", {
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

# ---- derive_var_atoxgr, test 12: CTCAEv5 Alanine aminotransferase increased ----
test_that("derive_var_atoxgr, test 12: CTCAEv5 Alanine aminotransferase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_alt_ctcv5_norm <- expected_alt_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      FLAG = "NORMAL"
    )

  # create records with baseline abnormal and apply criteria
  expected_alt_ctcv5_abn <- tibble::tribble(
    ~ATOXDSCH, ~AVAL, ~BASE, ~AVALU, ~ATOXGRH,
    "Not a term", 80, 40, NA_character_, NA,
    NA_character_, 60, 40, NA_character_, NA,
    "Alanine aminotransferase Increased", 801, 40, NA_character_, "4",
    "Alanine aminotransferase Increased", 800, 40, NA_character_, "3",
    "Alanine aminotransferase Increased", 201, 40, NA_character_, "3",
    "Alanine aminotransferase Increased", 200, 40, NA_character_, "2",
    "Alanine aminotransferase Increased", 121, 40, NA_character_, "2",
    "Alanine aminotransferase Increased", 120, 40, NA_character_, "1",
    "Alanine aminotransferase Increased", 60, 40, NA_character_, "1",
    "Alanine aminotransferase Increased", 59, 40, NA_character_, "0",
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased", 100, NA, NA_character_, NA,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased", NA, 40, NA_character_, NA,
  ) %>%
    # set BASE to be abnormal and create FLAG
    mutate(
      ANRHI = BASE - 1,
      FLAG = "ABNORMAL"
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
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alt_ctcv5,
    compare = actual_alt,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "FLAG", "AVALU")
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


# ---- derive_var_atoxgr, test 13: CTCAEv4 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr, test 13: CTCAEv4 Alkaline phosphatase increased", {
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

# ---- derive_var_atoxgr, test 14: CTCAEv5 Alkaline phosphatase increased ----
test_that("derive_var_atoxgr, test 14: CTCAEv5 Alkaline phosphatase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_alkph_ctcv5_norm <- expected_alkph_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      FLAG = "NORMAL"
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
    # set BASE to be abnormal and create FLAG
    mutate(
      ANRHI = BASE - 1,
      FLAG = "ABNORMAL"
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
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_alkph_ctcv5,
    compare = actual_alkph,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "FLAG", "AVALU")
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

# ---- derive_var_atoxgr, test 15: CTCAEv4 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr, test 15: CTCAEv4 Aspartate aminotransferase increased", {
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

# ---- derive_var_atoxgr, test 16: CTCAEv5 Aspartate aminotransferase increased ----
test_that("derive_var_atoxgr, test 16: CTCAEv5 Aspartate aminotransferase increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_ast_ctcv5_norm <- expected_ast_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      FLAG = "NORMAL"
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
    # set BASE to be abnormal and create FLAG
    mutate(
      ANRHI = BASE - 1,
      FLAG = "ABNORMAL"
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
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ast_ctcv5,
    compare = actual_ast,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "FLAG", "AVALU")
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

# ---- derive_var_atoxgr, test 17: CTCAEv4 Blood bilirubin increased ----
test_that("derive_var_atoxgr, test 17: CTCAEv4 Blood bilirubin increased", {
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

### Blood bilirubin increased
### NCICTCAEv5 same criteria as NCICTCAEv4 when BASELINE is normal
### Grade 4: >10.0 x ULN if BL was normal OR >10.0 x BL if BL was abnormal
### Grade 3: >3.0 - 10.0 x ULN if BL was normal OR >3.0 - 10.0 x BL
### Grade 2: >1.5 - 3.0 x ULN if BL was normal OR >1.5 - 3.0 x BL
### Grade 1: >ULN - 1.5 x ULN if BL was normal OR >1.0 - 1.5 x BL

# ---- derive_var_atoxgr, test 18: CTCAEv5  Blood bilirubin increased ----
test_that("derive_var_atoxgr, test 18: CTCAEv5  Blood bilirubin increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_bili_ctcv5_norm <- expected_bili_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      FLAG = "NORMAL"
    )

  # create records with abnormal BASE then add records with normal BASE
  expected_bili_ctcv5 <- expected_bili_ctcv4 %>%
    # set BASE to ANRHI then make ANRHI < BASE
    mutate(
      BASE = ANRHI,
      ANRHI = ANRHI - 1,
      FLAG = "ABNORMAL"
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
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_bili_ctcv5,
    compare = actual_bili,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "FLAG", "AVALU")
  )
})



### CD4 Lymphocytes decreased
### NCICTCAEv5 same criteria as NCICTCAEv4
### Grade 4: <0.05 x 10e9 /L
### Grade 3: <0.2 x 0.05 - 10e9 /L
### Grade 2: <0.5 - 0.2 x 10e9 /L
### Grade 1: <LLN - 0.5 x 10e9 /L

expected_cd4 <- tibble::tribble(
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

input_cd4 <- expected_cd4 %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 19: CTCAEv4 CD4 Lymphocytes decreased ----
test_that("derive_var_atoxgr, test 19: CTCAEv4 CD4 Lymphocytes decreased", {
  actual_cd4 <- derive_var_atoxgr_dir(
    input_cd4,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4,
    compare = actual_cd4,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 20: CTCAEv5 CD4 Lymphocytes decreased ----
test_that("derive_var_atoxgr, test 20: CTCAEv5 CD4 Lymphocytes decreased", {
  actual_cd4 <- derive_var_atoxgr_dir(
    input_cd4,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_cd4,
    compare = actual_cd4,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

### Cholesterol high
### NCICTCAEv5 same criteria as NCICTCAEv4
### Grade 4: >12.92 mmol/L
### Grade 3: >10.34 - 12.92 mmol/L
### Grade 2: >7.75 - 10.34 mmol/L
### Grade 1: >ULN - 7.75 mmol/L

expected_choles <- tibble::tribble(
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
input_choles <- expected_choles %>%
  select(-ATOXGRH)

# ---- derive_var_atoxgr, test 21: CTCAEv4 Cholesterol high ----
test_that("derive_var_atoxgr, test 21: CTCAEv4 Cholesterol high", {
  actual_choles <- derive_var_atoxgr_dir(
    input_choles,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles,
    compare = actual_choles,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 22: CTCAEv5 Cholesterol high ----
test_that("derive_var_atoxgr, test 22: CTCAEv5 Cholesterol high", {
  actual_choles <- derive_var_atoxgr_dir(
    input_choles,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_choles,
    compare = actual_choles,
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

# ---- derive_var_atoxgr, test 23: CTCAEv4 CPK increased ----
test_that("derive_var_atoxgr, test 23: CTCAEv4 CPK increased", {
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

# ---- derive_var_atoxgr, test 24: CTCAEv5 CPK increased ----
test_that("derive_var_atoxgr, test 24: CTCAEv5 CPK increased", {
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

# ---- derive_var_atoxgr, test 25: CTCAEv4 Creatinine increased ----
test_that("derive_var_atoxgr, test 25: CTCAEv4 Creatinine increased", {
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

# ---- derive_var_atoxgr, test 26: CTCAEv4 Creatinine increased ----
test_that("derive_var_atoxgr, test 26: CTCAEv4 Creatinine increased", {
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

### Fibrinogen decreased (NCICTCv4)
### Grade 4: <0.25 x LLN or 75% decrease from baseline or absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN or 50 - <75% decrease from baseline
### Grade 2: <0.75 - 0.5 x LLN or 25 - <50% decrease from baseline
### Grade 1: <1.0 - 0.75 x LLN or <25% decrease from baseline

# ---- derive_var_atoxgr, test 27: CTCAEv4 Fibrinogen decreased ----
test_that("derive_var_atoxgr, test 27: CTCAEv4 Fibrinogen decreased", {
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

### Fibrinogen decreased (NCICTCv5)
### Grade 4: <0.25 x LLN OR if abnormal, 75% dec. from BL OR absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN OR if abnormal, 50 - <75% dec. from BL
### Grade 2: <0.75 - 0.5 x LLN OR if abnormal, 25 - <50% dec. from BL
### Grade 1: <1.0 - 0.75 x LLN OR if abnormal, <25% dec. from BL

# ---- derive_var_atoxgr, test 28: CTCAEv5 Fibrinogen decreased ----
test_that("derive_var_atoxgr, test 28: CTCAEv5 Fibrinogen decreased", {
  expected_fib <- tibble::tribble(
    ~ATOXDSCL, ~AVAL, ~ANRLO, ~PCHG, ~AVALU, ~ATOXGRL,
    "Not a term", 9, 10, 40, "g/L", NA,
    NA_character_, 10, 10, 40, "g/L", NA,
    # Satisfies < 0.5 for grade 4 - other criteria missing
    "Fibrinogen decreased", 0.49, NA, NA, "g/L", "4",
    # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
    "Fibrinogen decreased", 0.49, 1, -51, "g/L", "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG missing
    "Fibrinogen decreased", 0.5, 2.1, NA, "g/L", "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
    "Fibrinogen decreased", 0.5, 2.1, -51, "g/L", "4",
    # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
    "Fibrinogen decreased", 0.5, 1.1, -75, "g/L", "4",
    # Satisfies < 0.5*LLN for grade 3 - PCHG missing
    "Fibrinogen decreased", 1, 2.1, NA, "g/L", "3",
    # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
    "Fibrinogen decreased", 1, 2.1, -49, "g/L", "3",
    # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
    "Fibrinogen decreased", 1, 2, -50, "g/L", "3",
    # Satisfies < 0.75*LLN for grade 2 - PCHG missing
    "Fibrinogen decreased", 1.5, 2.1, NA, "g/L", "2",
    # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
    "Fibrinogen decreased", 1.5, 2.1, -10, "g/L", "2",
    # Satisfies <=25% for grade 2 - LLN satisfies grade 1
    "Fibrinogen decreased", 1.5, 1.6, -25, "g/L", "2",
    # Satisfies < LLN for grade 1 - PCHG missing
    "Fibrinogen decreased", 2, 2.1, NA, "g/L", "1",
    # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
    "Fibrinogen decreased", 2, 2.1, 10, "g/L", "1",
    # Satisfies grade 0 - AVAL >= LLN AND no % descrease
    "Fibrinogen decreased", 1.5, 1.5, 0, "g/L", "0",
    # Satisfies % decrease for grade 1 - AVAL = LLN so not abnormal
    "Fibrinogen decreased", 1.5, 1.5, -1, "g/L", "0",
    # AVAL >= LLN - PCT missing but its normal so ignore PCT
    "Fibrinogen decreased", 1.5, 1.5, NA, "g/L", "0",
    # Satisfies <=75% decrease for grade 4 - LLN missing do not know its abnormal
    "Fibrinogen decreased", 1, NA, -75, "g/L", NA,
    # Satisfies <=50% decrease for grade 3 - LLN missing do not know its abnormal
    "Fibrinogen decreased", 1, NA, -50, "g/L", NA,
    # Satisfies <=25% decrease for grade 2 - LLN missing do not know its abnormal
    "Fibrinogen decreased", 1.5, NA, -25, "g/L", NA,
    # Satisfies % decrease for grade 1 - LLN missing do not know its abnormal
    "Fibrinogen decreased", 1.5, NA, -1, "g/L", NA,
    # PCT >= 0 BUT LLN missing cannot grade as NORMAL
    "Fibrinogen decreased", 1.5, NA, 10, "g/L", NA,
    # AVAL missing cannot grade
    "Fibrinogen decreased", NA, 1.5, 10, "g/L", NA,
    # wrong unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased", 1.5, 1.5, 0, "g/dL", NA,
    # missing unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased", 1.5, 1.5, 0, NA, NA,
  )
  input_fib <- expected_fib %>%
    select(-ATOXGRL)

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

# ---- derive_var_atoxgr, test 29: CTCAEv4 GGT increased ----
test_that("derive_var_atoxgr, test 29: CTCAEv4 GGT increased", {
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

# ---- derive_var_atoxgr, test 30: CTCAEv5 GGT increased ----
test_that("derive_var_atoxgr, test 30: CTCAEv5 GGT increased", {
  # V5 and V4 criteria identical when BASELINE normal
  expected_ggt_ctcv5_norm <- expected_ggt_ctcv4 %>%
    # set BASE to be normal and create FLAG
    mutate(
      BASE = ANRHI,
      FLAG = "NORMAL"
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
    # set BASE to be abnormal and create FLAG
    mutate(
      ANRHI = BASE - 1,
      FLAG = "ABNORMAL"
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
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_ggt_ctcv5,
    compare = actual_ggt,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "FLAG", "AVALU")
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

# ---- derive_var_atoxgr, test 31: CTCAEv4 Haptoglobin decreased ----
test_that("derive_var_atoxgr, test 31: CTCAEv4 Haptoglobin decreased", {
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

# ---- derive_var_atoxgr, test 32: CTCAEv5 Haptoglobin decreased ----
test_that("derive_var_atoxgr, test 32: CTCAEv5 Haptoglobin decreased", {
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

### Hemoglobin increased
# NCICTCAEv5 same as NCICTCAEv4 when BASE is normal
### Grade 3: Increase in >4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 2: Increase in >2 - 4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 1: Increase in >0 - 2 gm/dL above ULN or above baseline if baseline is above ULN

expected_hgbi <- tibble::tribble(
  ~ATOXDSCH, ~AVAL, ~BASE, ~ANRHI, ~AVALU, ~ATOXGRH, ~TESTNUM, ~V5,
  "Not a term", 80, 120, 200, "g/L", NA, 1, "Y",
  NA_character_, 60, 50, 100, "g/L", NA, 2, "Y",
  # BASE greater than ANRHI
  "Hemoglobin increased", 106, 65, 60, "g/L", "3", 3, "N",
  "Hemoglobin increased", 105, 65, 60, "g/L", "2", 4, "N",
  "Hemoglobin increased", 86, 65, 60, "g/L", "2", 5, "N",
  "Hemoglobin increased", 85, 65, 60, "g/L", "1", 6, "N",
  "Hemoglobin increased", 66, 65, 60, "g/L", "1", 7, "N",
  "Hemoglobin increased", 65, 65, 60, "g/L", "0", 8, "N",
  "Hemoglobin increased", NA, 65, 60, "g/L", NA, 9, "N",
  # BASE less than or equal to ANRHI
  "Hemoglobin increased", 106, 60, 65, "g/L", "3", 10, "Y",
  "Hemoglobin increased", 105, 60, 65, "g/L", "2", 11, "Y",
  "Hemoglobin increased", 86, 60, 65, "g/L", "2", 12, "Y",
  "Hemoglobin increased", 85, 60, 65, "g/L", "1", 13, "Y",
  "Hemoglobin increased", 66, 60, 65, "g/L", "1", 14, "Y",
  "Hemoglobin increased", 65, 60, 65, "g/L", "0", 15, "Y",
  "Hemoglobin increased", NA, 60, 65, "g/L", NA, 16, "Y",
  # BASE missing
  "Hemoglobin increased", 106, NA, 65, "g/L", "3", 17, "N",
  "Hemoglobin increased", 105, NA, 65, "g/L", "2", 18, "N",
  "Hemoglobin increased", 86, NA, 65, "g/L", "2", 19, "N",
  "Hemoglobin increased", 85, NA, 65, "g/L", "1", 20, "N",
  "Hemoglobin increased", 66, NA, 65, "g/L", "1", 21, "N",
  "Hemoglobin increased", 65, NA, 65, "g/L", "0", 22, "N",
  "Hemoglobin increased", NA, NA, 65, "g/L", NA, 23, "N",
  # Unit missing cannot grade
  "Hemoglobin increased", 200, 61, 65, NA, NA, 24, "Y",
  # ANRHI missing - cannot grade
  "Hemoglobin increased", 200, 60, NA, "g/L", NA, 25, "Y",
  # AVAL missing cannot grade
  "Hemoglobin increased", NA, 60, 65, "g/L", NA, 26, "Y",
)

# ---- derive_var_atoxgr, test 33: CTCAEv4 Hemoglobin increased ----
test_that("derive_var_atoxgr, test 33: CTCAEv4 Hemoglobin increased", {
  input_hgbi <- expected_hgbi %>%
    select(-ATOXGRH)

  actual_hgbi <- derive_var_atoxgr_dir(
    input_hgbi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi,
    compare = actual_hgbi,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 34: CTCAEv5 Hemoglobin increased ----
test_that("derive_var_atoxgr, test 34: CTCAEv5 Hemoglobin increased", {
  expected_hgbi <- expected_hgbi %>%
    filter(V5 == "Y")

  input_hgbi <- expected_hgbi %>%
    select(-ATOXGRH)

  actual_hgbi <- derive_var_atoxgr_dir(
    input_hgbi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_hgbi,
    compare = actual_hgbi,
    keys = c("TESTNUM")
  )
})

### INR increased (NCICTCV4)
### NCICTCV5 different for grade 1
### Grade 3: >2.5 x ULN; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5 x ULN; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1 - 1.5 x ULN; >1 - 1.5 times above baseline if on anticoagulation

# ---- derive_var_atoxgr, test 35: CTCAEv4 INR increased ----
test_that("derive_var_atoxgr, test 35: CTCAEv4 INR increased", {
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

### INR increased (NCICTCV4)
### NCICTCV5 different for grade 1
### Grade 3: >2.5 x ULN; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5 x ULN; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1.2 - 1.5 x ULN; >1 - 1.5 times above baseline if on anticoagulation

# ---- derive_var_atoxgr, test 36: CTCAEv5 INR increased ----
test_that("derive_var_atoxgr, test 36: CTCAEv5 INR increased", {
  expected_inri <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~BASE, ~ANRHI, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    120,   200,    NA_character_, NA,       1,
    NA_character_,   60,    50,    100,    NA_character_, NA,       2,
    # GRADE derived from AVAL against ANRHI
    "INR IncreaSed", 251,   200,   100,    NA_character_, "3",      3,
    "INR Increased", 250,   199,   100,    NA_character_, "2",      4,
    "INR Increased", 151,   150,   100,    NA_character_, "2",      5,
    "INR Increased", 150,   150,   100,    NA_character_, "1",      6,
    "INR Increased", 121,   150,   100,    NA_character_, "1",      7,
    "INR Increased", 120,   120,   100,    NA_character_, "0",      8,
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

# ---- derive_var_atoxgr, test 37: CTCAEv4 Lipase increased ----
test_that("derive_var_atoxgr, test 37: CTCAEv4 Lipase increased", {
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

# ---- derive_var_atoxgr, test 38: CTCAEv5 Lipase increased ----
test_that("derive_var_atoxgr, test 38: CTCAEv5 Lipase increased", {
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


### Lymphocyte count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <0.2 x 10e9 /L
### Grade 3: <0.5 - 0.2 x 10e9 /L
### Grade 2: <0.8 - 0.5 x 10e9 /L
### Grade 1: <LLN - 0.8 x 10e9/L

expected_lymd <- tibble::tribble(
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
input_lymd <- expected_lymd %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 39: CTCAEv4 Lymphocyte count decreased ----
test_that("derive_var_atoxgr, test 39: CTCAEv4 Lymphocyte count decreased", {
  actual_lymd <- derive_var_atoxgr_dir(
    input_lymd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd,
    compare = actual_lymd,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 40: CTCAEv5 Lymphocyte count decreased ----
test_that("derive_var_atoxgr, test 40: CTCAEv5 Lymphocyte count decreased", {
  actual_lymd <- derive_var_atoxgr_dir(
    input_lymd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymd,
    compare = actual_lymd,
    keys = c("TESTNUM")
  )
})

### Lymphocyte count increased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: >20,000/mm3
### Grade 2: >4000/mm3 - 20,000/mm3

expected_lymi <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,
  "Lymphocyte count increased", 21,     NA,     NA,     "10^9/L",  "3",
  "Lymphocyte count increased", 20,     NA,     NA,     "10^9/L",  "2",
  "Lymphocyte count increased", 4.1,    NA,     NA,     "10^9/L",  "2",
  "Lymphocyte count increased", 4,      NA,     NA,     "10^9/L",  "0",
  # Unit missing cannot grade
  "Lymphocyte count increased", 4,      NA,     NA,     NA,        NA,
  # AVAL missing cannot grade
  "Lymphocyte count increased", NA,     NA,     NA,     "10^9/L",  NA,
)
input_lymi <- expected_lymi %>%
  select(-ATOXGRH)

# ---- derive_var_atoxgr, test 41: CTCAEv4 Lymphocyte count increased ----
test_that("derive_var_atoxgr, test 41: CTCAEv4 Lymphocyte count increased", {
  actual_lymi <- derive_var_atoxgr_dir(
    input_lymi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi,
    compare = actual_lymi,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 42: CTCAEv5 Lymphocyte count increased ----
test_that("derive_var_atoxgr, test 42: CTCAEv5 Lymphocyte count increased", {
  actual_lymi <- derive_var_atoxgr_dir(
    input_lymi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_lymi,
    compare = actual_lymi,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Neutrophil count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <1.0 - 0.5 x 10e9 /L
### Grade 2: <1.5 - 1.0 x 10e9 /L
### Grade 1: <LLN - 1.5 x 10e9 /L

expected_neut <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,
  # ANRLO not missing
  "Neutrophil count decreased", 0.49,   2,      NA,     "10^9/L",  "4",
  "Neutrophil count decreased", 0.5,    2,      NA,     "10^9/L",  "3",
  "Neutrophil count decreased", 0.99,   2,      NA,     "10^9/L",  "3",
  "Neutrophil count decreased", 1,      2,      NA,     "10^9/L",  "2",
  "Neutrophil count decreased", 1.49,   2,      NA,     "10^9/L",  "2",
  "Neutrophil count decreased", 1.5,    2,      NA,     "10^9/L",  "1",
  "Neutrophil count decreased", 1.9,    2,      NA,     "10^9/L",  "1",
  "Neutrophil count decreased", 2,      2,      NA,     "10^9/L",  "0",
  # ANRLO missing - can grade 2-4
  "Neutrophil count decreased", 0.49,   NA,     NA,     "10^9/L",  "4",
  "Neutrophil count decreased", 0.5,    NA,     NA,     "10^9/L",  "3",
  "Neutrophil count decreased", 0.99,   NA,     NA,     "10^9/L",  "3",
  "Neutrophil count decreased", 1,      NA,     NA,     "10^9/L",  "2",
  "Neutrophil count decreased", 1.49,   NA,     NA,     "10^9/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Neutrophil count decreased", 1.5,    NA,     NA,     "10^9/L",  NA,
  "Neutrophil count decreased", 1.9,    NA,     NA,     "10^9/L",  NA,
  "Neutrophil count decreased", 2,      NA,     NA,     "10^9/L",  NA,
  # Unit missing cannot grade
  "Neutrophil count decreased", 2,      2,      NA,     NA,        NA,
  # AVAL missing cannot grade
  "Neutrophil count decreased", NA,     2,      NA,     "10^9/L",  NA,
)
input_neut <- expected_neut %>%
  select(-ATOXGRL)


# ---- derive_var_atoxgr, test 43: CTCAEv4 Neutrophil count decreased ----
test_that("derive_var_atoxgr, test 43: CTCAEv4 Neutrophil count decreased", {
  actual_neut <- derive_var_atoxgr_dir(
    input_neut,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut,
    compare = actual_neut,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 44: CTCAEv5 Neutrophil count decreased ----
test_that("derive_var_atoxgr, test 44: CTCAEv5 Neutrophil count decreased", {
  actual_neut <- derive_var_atoxgr_dir(
    input_neut,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_neut,
    compare = actual_neut,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Platelet count decreased
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <50.0 - 25.0 x 10e9 /L
### Grade 2: <75.0 - 50.0 x 10e9 /L
### Grade 1: <LLN - 75.0 x 10e9 /L

expected_plate <- tibble::tribble(
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
input_plate <- expected_plate %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 45: CTCAEv4 Platelet count decreased ----
test_that("derive_var_atoxgr, test 45: CTCAEv4 Platelet count decreased", {
  actual_plate <- derive_var_atoxgr_dir(
    input_plate,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate,
    compare = actual_plate,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 46: CTCAEv5 Platelet count decreased ----
test_that("derive_var_atoxgr, test 46: CTCAEv5 Platelet count decreased", {
  actual_plate <- derive_var_atoxgr_dir(
    input_plate,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_plate,
    compare = actual_plate,
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

# ---- derive_var_atoxgr, test 47: CTCAEv4 Serum amylase increased ----
test_that("derive_var_atoxgr, test 47: CTCAEv4 Serum amylase increased", {
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

# ---- derive_var_atoxgr, test 48: CTCAEv5 Serum amylase increased ----
test_that("derive_var_atoxgr, test 48: CTCAEv5 Serum amylase increased", {
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

# ---- derive_var_atoxgr, test 49: CTCAEv4 White blood cell decreased ----
test_that("derive_var_atoxgr, test 49: CTCAEv4 White blood cell decreased", {
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

# ---- derive_var_atoxgr, test 50: CTCAEv5 White blood cell decreased ----
test_that("derive_var_atoxgr, test 50: CTCAEv5 White blood cell decreased", {
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


## Metabolism and nutrition disorders

### Hypercalcemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >3.4 mmol/L
### Grade 3: >3.1 - 3.4 mmol/L
### Grade 2: >2.9 - 3.1 mmol/L
### Grade 1: >ULN - 2.9 mmol/L

expected_calci <- tibble::tribble(
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
input_calci <- expected_calci %>%
  select(-ATOXGRH)

# ---- derive_var_atoxgr, test 51: CTCAEv4 Hypercalcemia ----
test_that("derive_var_atoxgr, test 51: CTCAEv4 Hypercalcemia", {
  actual_calci <- derive_var_atoxgr_dir(
    input_calci,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci,
    compare = actual_calci,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 52: CTCAEv5 Hypercalcemia ----
test_that("derive_var_atoxgr, test 52: CTCAEv5 Hypercalcemia", {
  actual_calci <- derive_var_atoxgr_dir(
    input_calci,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calci,
    compare = actual_calci,
    keys = c("TESTNUM")
  )
})

### Hypercalcemia (Ionized)
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >1.8 mmol/L
### Grade 3: >1.6 - 1.8 mmol/L
### Grade 2: >1.5 - 1.6 mmol/L
### Grade 1: >ULN - 1.5 mmol/L

expected_calioni <- tibble::tribble(
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
input_calioni <- expected_calioni %>%
  select(-ATOXGRH)


# ---- derive_var_atoxgr, test 53: CTCAEv4 Hypercalcemia (Ionized) ----
test_that("derive_var_atoxgr, test 53: CTCAEv4 Hypercalcemia (Ionized)", {
  actual_calioni <- derive_var_atoxgr_dir(
    input_calioni,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni,
    compare = actual_calioni,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 54: CTCAEv5 Hypercalcemia (Ionized) ----
test_that("derive_var_atoxgr, test 54: CTCAEv5 Hypercalcemia (Ionized)", {
  actual_calioni <- derive_var_atoxgr_dir(
    input_calioni,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calioni,
    compare = actual_calioni,
    keys = c("TESTNUM")
  )
})

### Hyperglycemia (Fasting) (NCICTCAEv4)
### not included in NCICTCAEv5
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L
### Grade 2: >8.9 - 13.9 mmol/L
### Grade 1: >ULN - 8.9 mmol/L

# ---- derive_var_atoxgr, test 55: CTCAEv4 Hyperglycemia (Fasting) ----
test_that("derive_var_atoxgr, test 55: CTCAEv4 Hyperglycemia (Fasting)", {
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

### Hyperglycemia (NCICTCAEv4)
### not included in NCICTCAEv5
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L

# ---- derive_var_atoxgr, test 56: CTCAEv4 Hyperglycemia ----
test_that("derive_var_atoxgr, test 56: CTCAEv4 Hyperglycemia", {
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

### Hyperkalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
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

# ---- derive_var_atoxgr, test 57: CTCAEv4 Hyperkalemia ----
test_that("derive_var_atoxgr, test 57: CTCAEv4 Hyperkalemia", {
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

# ---- derive_var_atoxgr, test 58: CTCAEv5 Hyperkalemia ----
test_that("derive_var_atoxgr, test 58: CTCAEv5 Hyperkalemia", {
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

### Hypermagnesemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >3.30 mmol/L
### Grade 3: >1.23 - 3.30 mmol/L
### Grade 1: >ULN - 1.23 mmol/L

expected_magni <- tibble::tribble(
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
input_magni <- expected_magni %>%
  select(-ATOXGRH)


# ---- derive_var_atoxgr, test 59: CTCAEv4 Hypermagnesemia ----
test_that("derive_var_atoxgr, test 59: CTCAEv4 Hypermagnesemia", {
  actual_magni <- derive_var_atoxgr_dir(
    input_magni,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni,
    compare = actual_magni,
    keys = c("TESTNUM")
  )
})

# ---- derive_var_atoxgr, test 60: CTCAEv5 Hypermagnesemia ----
test_that("derive_var_atoxgr, test 60: CTCAEv5 Hypermagnesemia", {
  actual_magni <- derive_var_atoxgr_dir(
    input_magni,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magni,
    compare = actual_magni,
    keys = c("TESTNUM")
  )
})

### Hypernatremia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
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

# ---- derive_var_atoxgr, test 61: CTCAEv4 Hypernatremia ----
test_that("derive_var_atoxgr, test 61: CTCAEv4 Hypernatremia", {
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

# ---- derive_var_atoxgr, test 62: CTCAEv5 Hypernatremia ----
test_that("derive_var_atoxgr, test 62: CTCAEv5 Hypernatremia", {
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

### Hypertriglyceridemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: >11.4 mmol/L
### Grade 3: >5.7 mmol/L - 11.4 mmol/L
### Grade 2: >3.42 mmol/L - 5.7 mmol/L
### Grade 1: 1.71 mmol/L - 3.42 mmol/L

expected_trigi <- tibble::tribble(
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
input_trigi <- expected_trigi %>%
  select(-ATOXGRH)

# ---- derive_var_atoxgr, test 63: CTCAEv4 Hypertriglyceridemia ----
test_that("derive_var_atoxgr, test 63: CTCAEv4 Hypertriglyceridemia", {
  actual_trigi <- derive_var_atoxgr_dir(
    input_trigi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi,
    compare = actual_trigi,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 64: CTCAEv5 Hypertriglyceridemia ----
test_that("derive_var_atoxgr, test 64: CTCAEv5 Hypertriglyceridemia", {
  actual_trigi <- derive_var_atoxgr_dir(
    input_trigi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_trigi,
    compare = actual_trigi,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

### Hyperuricemia (NCICTCAEv4)
### NCICTCAEv5 only has grade 3
### Grade 4: >0.59 mmol/L;
### Grade 3: >ULN - 10 mg/dL (0.59 mmol/L)

expected_albi <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
  "Not a term",     591,    0,      200,    "umol/L",  NA,
  NA_character_,    591,    0,      200,    "umol/L",  NA,
  # ANRHI not missing
  "Hyperuricemia",  591,    0,      200,    "umol/L",  "4",
  "Hyperuricemia",  590,    0,      200,    "umol/L",  "3",
  "Hyperuricemia",  201,    0,      200,    "umol/L",  "3",
  "Hyperuricemia",  200,    0,      200,    "umol/L",  "0",
  # ANRHI missing - can grade 4
  "Hyperuricemia",  591,    0,      NA,     "umol/L",  "4",
  # ANRHI missing - can NOT grade 0 or 3
  "Hyperuricemia",  590,    0,      NA,     "umol/L",  NA,
  "Hyperuricemia",  201,    0,      NA,     "umol/L",  NA,
  "Hyperuricemia",  200,    0,      NA,     "umol/L",  NA,
  # Unit missing cannot grade
  "Hyperuricemia",  200,    0,      200,    NA,        NA,
  # AVAL missing cannot grade
  "Hyperuricemia",  NA,     0,      200,    "umol/L",  NA,
)
input_albi <- expected_albi %>%
  select(-ATOXGRH)

### Hyperuricemia (NCICTCAEv5)
### NCICTCAEv5 only has grade 3
### Grade 3: >ULN

# ---- derive_var_atoxgr, test 65: CTCAEv5 Hyperuricemia ----
test_that("derive_var_atoxgr, test 65: CTCAEv5 Hyperuricemia", {
  expected_albi <- expected_albi %>%
    filter(ATOXGRH != "4")
  input_albi <- expected_albi %>%
    select(-ATOXGRH)

  actual_albi <- derive_var_atoxgr_dir(
    input_albi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albi,
    compare = actual_albi,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 66: CTCAEv4 Hyperuricemia ----
test_that("derive_var_atoxgr, test 66: CTCAEv4 Hyperuricemia", {
  actual_albi <- derive_var_atoxgr_dir(
    input_albi,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albi,
    compare = actual_albi,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### Hypoalbuminemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 3: <20 g/L
### Grade 2: <30 - 20 g/L
### Grade 1: <LLN - 30 g/L

expected_albd <- tibble::tribble(
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
input_albd <- expected_albd %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 67: CTCAEv4 Hypoalbuminemia ----
test_that("derive_var_atoxgr, test 67: CTCAEv4 Hypoalbuminemia", {
  actual_albd <- derive_var_atoxgr_dir(
    input_albd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd,
    compare = actual_albd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 68: CTCAEv5 Hypoalbuminemia ----
test_that("derive_var_atoxgr, test 68: CTCAEv5 Hypoalbuminemia", {
  actual_albd <- derive_var_atoxgr_dir(
    input_albd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_albd,
    compare = actual_albd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hypocalcemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <1.5 mmol/L
### Grade 3: <1.75 - 1.5 mmol/L
### Grade 2: <2.0 - 1.75 mmol/L
### Grade 1: <LLN - 2.0 mmol/L

expected_calcd <- tibble::tribble(
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
input_calcd <- expected_calcd %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 69: CTCAEv4 Hypocalcemia ----
test_that("derive_var_atoxgr, test 69: CTCAEv4 Hypocalcemia", {
  actual_calcd <- derive_var_atoxgr_dir(
    input_calcd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd,
    compare = actual_calcd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 70: CTCAEv5 Hypocalcemia ----
test_that("derive_var_atoxgr, test 70: CTCAEv5 Hypocalcemia", {
  actual_calcd <- derive_var_atoxgr_dir(
    input_calcd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_calcd,
    compare = actual_calcd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hypocalcemia (Ionized)
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <0.8 mmol/L
### Grade 3: <0.9 - 0.8 mmol/L
### Grade 2: <1.0 - 0.9 mmol/L
### Grade 1: <LLN - 1.0 mmol/L

expected_caliond <- tibble::tribble(
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
input_caliond <- expected_caliond %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 71: CTCAEv4 Hypocalcemia (Ionized) ----
test_that("derive_var_atoxgr, test 71: CTCAEv4 Hypocalcemia (Ionized)", {
  actual_caliond <- derive_var_atoxgr_dir(
    input_caliond,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond,
    compare = actual_caliond,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 72: CTCAEv5 Hypocalcemia (Ionized) ----
test_that("derive_var_atoxgr, test 72: CTCAEv5 Hypocalcemia (Ionized)", {
  actual_caliond <- derive_var_atoxgr_dir(
    input_caliond,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_caliond,
    compare = actual_caliond,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hypoglycemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <1.7 mmol/L
### Grade 3: <2.2 - 1.7 mmol/L
### Grade 2: <3.0 - 2.2 mmol/L
### Grade 1: <LLN - 3.0 mmol/L

expected_glycd <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
  "Not a term",    1.69,   4,      100,    "mmol/L",  NA,
  NA_character_,   1.69,   4,      100,    "mmol/L",  NA,
  # ANRLO not missing
  "Hypoglycemia",  1.69,   4,      100,    "mmol/L",  "4",
  "Hypoglycemia",  1.7,    4,      100,    "mmol/L",  "3",
  "Hypoglycemia",  2.19,   4,      100,    "mmol/L",  "3",
  "Hypoglycemia",  2.2,    4,      100,    "mmol/L",  "2",
  "Hypoglycemia",  2.9,    4,      100,    "mmol/L",  "2",
  "Hypoglycemia",  3,      4,      100,    "mmol/L",  "1",
  "Hypoglycemia",  3.9,    4,      100,    "mmol/L",  "1",
  "Hypoglycemia",  4,      4,      100,    "mmol/L",  "0",
  # ANRLO missing - can grade 2-4
  "Hypoglycemia",  1.69,   NA,     100,    "mmol/L",  "4",
  "Hypoglycemia",  1.7,    NA,     100,    "mmol/L",  "3",
  "Hypoglycemia",  2.19,   NA,     100,    "mmol/L",  "3",
  "Hypoglycemia",  2.2,    NA,     100,    "mmol/L",  "2",
  "Hypoglycemia",  2.9,    NA,     100,    "mmol/L",  "2",
  # ANRLO missing - can NOT grade 0 or 1
  "Hypoglycemia",  3,      NA,     100,    "mmol/L",  NA,
  "Hypoglycemia",  3.9,    NA,     100,    "mmol/L",  NA,
  "Hypoglycemia",  4,      NA,     100,    "mmol/L",  NA,
  # Unit missing cannot grade
  "Hypoglycemia",  4,      4,      100,    NA,        NA,
  # AVAL missing cannot grade
  "Hypoglycemia",  NA,     4,      100,    "mmol/L",  NA,
)
input_glycd <- expected_glycd %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 73: CTCAEv4 Hypoglycemia ----
test_that("derive_var_atoxgr, test 73: CTCAEv4 Hypoglycemia", {
  actual_glycd <- derive_var_atoxgr_dir(
    input_glycd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd,
    compare = actual_glycd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 74: CTCAEv5 Hypoglycemia ----
test_that("derive_var_atoxgr, test 74: CTCAEv5 Hypoglycemia", {
  actual_glycd <- derive_var_atoxgr_dir(
    input_glycd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_glycd,
    compare = actual_glycd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hypokalemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
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

# ---- derive_var_atoxgr, test 75: CTCAEv4 Hypokalemia ----
test_that("derive_var_atoxgr, test 75: CTCAEv4 Hypokalemia", {
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

# ---- derive_var_atoxgr, test 76: CTCAEv5 Hypokalemia ----
test_that("derive_var_atoxgr, test 76: CTCAEv5 Hypokalemia", {
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

### Hypomagnesemia
### NCICTCAEv4 and NCICTCAEv5 criteria is the same
### Grade 4: <0.3 mmol/L
### Grade 3: <0.4 - 0.3 mmol/L
### Grade 2: <0.5 - 0.4 mmol/L
### Grade 1: <LLN - 0.5 mmol/L

expected_magnd <- tibble::tribble(
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
input_magnd <- expected_magnd %>%
  select(-ATOXGRL)

# ---- derive_var_atoxgr, test 77: CTCAEv4 Hypomagnesemia ----
test_that("derive_var_atoxgr, test 77: CTCAEv4 Hypomagnesemia", {
  actual_magnd <- derive_var_atoxgr_dir(
    input_magnd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd,
    compare = actual_magnd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

# ---- derive_var_atoxgr, test 78: CTCAEv5 Hypomagnesemia ----
test_that("derive_var_atoxgr, test 78: CTCAEv5 Hypomagnesemia", {
  actual_magnd <- derive_var_atoxgr_dir(
    input_magnd,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv5,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected_magnd,
    compare = actual_magnd,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### Hyponatremia (NCICTCAEv4)
### NCICTCAEv4 and NCICTCAEv5 essentially the same (slightly different text)
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

# ---- derive_var_atoxgr, test 79: CTCAEv4 Hyponatremia ----
test_that("derive_var_atoxgr, test 79: CTCAEv4 Hyponatremia", {
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

# ---- derive_var_atoxgr, test 80: CTCAEv5 Hyponatremia ----
test_that("derive_var_atoxgr, test 80: CTCAEv5 Hyponatremia", {
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

### Hypophosphatemia
### Only in NCICTCAEv4 (not NCICTCAEv5)
### Grade 4: <0.3 mmol/L
### Grade 3: <0.6 - 0.3 mmol/L
### Grade 2: <0.8 - 0.6 mmol/L
### Grade 1: <LLN - 0.8 mmol/L

# ---- derive_var_atoxgr, test 81: CTCAEv4 Hypophosphatemia ----
test_that("derive_var_atoxgr, test 81: CTCAEv4 Hypophosphatemia", {
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
