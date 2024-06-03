get_smq <- function(basket_select,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  if (basket_select$scope == "NARROW") {
    end <- 3
  } else {
    end <- 5
  }

  if (is.null(basket_select$name)) {
    basket_select$name <- paste("SMQ name of", basket_select$id)
  }
  terms <- tibble(TERMCHAR = paste(basket_select$name, "Term", c(1:end), "(", version, ")"))
  terms <- mutate(terms, SRCVAR = "AEDECOD", GRPNAME = basket_select$name)
  if (keep_id) {
    mutate(terms, GRPID = 42)
  } else {
    terms
  }
}

get_sdg <- function(basket_select,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  terms <- tibble(TERMCHAR = paste(basket_select$name, "Term", c(1:4)))
  terms <- mutate(terms, SRCVAR = "CMDECOD", GRPNAME = basket_select$name)
  if (keep_id) {
    mutate(terms, GRPID = 42)
  } else {
    terms
  }
}

cqterms <- tibble::tribble(
  ~TERMCHAR, ~TERMNUM,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L
) %>%
  mutate(SRCVAR = "AEDECOD")

# create_query_data ----
## Test 1: customized query defined by terms ----
test_that("create_query_data Test 1: customized query defined by terms", {
  cq <- query(
    prefix = "CQ01",
    name = "Application Site Issues",
    definition = cqterms
  )

  actual_output <- create_query_data(queries = list(cq))

  expected_output <- cqterms %>% mutate(
    GRPNAME = "Application Site Issues",
    PREFIX = "CQ01"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

## Test 2: customized query defined by SMQs ----
test_that("create_query_data Test 2: customized query defined by SMQs", {
  cq <- query(
    prefix = "CQ02",
    name = "Immune-Mediated Meningoencephalitis",
    definition = list(
      basket_select(
        name = "Noninfectious meningitis",
        scope = "NARROW",
        type = "smq"
      ),
      basket_select(
        name = "Noninfectious encephalitis",
        scope = "BROAD",
        type = "smq"
      )
    )
  )

  actual_output <- create_query_data(
    queries = list(cq),
    version = "20.0",
    get_terms_fun = get_smq
  )

  expected_output <-
    bind_rows(
      get_smq(
        basket_select(
          name = "Noninfectious meningitis",
          scope = "NARROW",
          type = "smq"
        ),
        version = "20.0"
      ),
      get_smq(
        basket_select(
          name = "Noninfectious encephalitis",
          scope = "BROAD",
          type = "smq"
        ),
        version = "20.0"
      )
    ) %>%
    mutate(
      GRPNAME = "Immune-Mediated Meningoencephalitis",
      PREFIX = "CQ02",
      VERSION = "20.0"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

## Test 3: customized query defined by terms and SMQs ----
test_that("create_query_data Test 3: customized query defined by terms and SMQs", {
  cq <- query(
    prefix = "CQ03",
    name = "Immune-Mediated Meningoencephalitis or Application Site Issues",
    definition = list(
      basket_select(
        name = "Noninfectious meningitis",
        scope = "NARROW",
        type = "smq"
      ),
      cqterms,
      basket_select(
        name = "Noninfectious encephalitis",
        scope = "BROAD",
        type = "smq"
      )
    )
  )

  actual_output <- create_query_data(
    queries = list(cq),
    version = "20.1",
    get_terms_fun = get_smq
  )

  expected_output <-
    bind_rows(
      get_smq(
        basket_select(
          name = "Noninfectious meningitis",
          scope = "NARROW",
          type = "smq"
        ),
        version = "20.1"
      ),
      cqterms,
      get_smq(
        basket_select(
          name = "Noninfectious encephalitis",
          scope = "BROAD",
          type = "smq"
        ),
        version = "20.1"
      )
    ) %>%
    mutate(
      GRPNAME = "Immune-Mediated Meningoencephalitis or Application Site Issues",
      PREFIX = "CQ03",
      VERSION = "20.1"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

## Test 4: SMQs ----
test_that("create_query_data Test 4: SMQs", {
  pregsmq <- query(
    prefix = "SMQ02",
    id = 13,
    add_scope_num = TRUE,
    definition = basket_select(
      name = "Pregnancy and neonatal topics (SMQ)",
      scope = "NARROW",
      type = "smq"
    )
  )

  pneuaegt <- query(
    prefix = "SMQ04",
    definition = basket_select(
      id = 8050L,
      scope = "BROAD",
      type = "smq"
    )
  )

  actual_output <-
    create_query_data(
      queries = list(pregsmq, pneuaegt),
      version = "20.0",
      get_terms_fun = get_smq
    )

  expected_output <-
    bind_rows(
      get_smq(
        basket_select(
          name = "Pregnancy and neonatal topics (SMQ)",
          scope = "NARROW",
          type = "smq"
        ),
        version = "20.0"
      ) %>%
        mutate(
          GRPNAME = "Pregnancy and neonatal topics (SMQ)",
          GRPID = 13,
          SCOPE = "NARROW",
          SCOPEN = 2,
          PREFIX = "SMQ02"
        ),
      get_smq(
        basket_select(
          id = 8050L,
          scope = "BROAD",
          type = "smq"
        ),
        version = "20.0"
      ) %>%
        mutate(
          SCOPE = "BROAD",
          PREFIX = "SMQ04"
        )
    ) %>%
    mutate(
      VERSION = "20.0"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

## Test 5: issues error if SMQs without meddra_version are requested ----
test_that("create_query_data Test 5: issues error if SMQs without meddra_version are requested", {
  pregsmq <- query(
    prefix = "SMQ02",
    definition = basket_select(
      name = "Pregnancy and neonatal topics (SMQ)",
      scope = "NARROW",
      type = "smq"
    )
  )

  expect_snapshot(
    create_query_data(
      queries = list(pregsmq),
      get_terms_fun = get_smq
    ),
    error = TRUE
  )
})

## Test 6: SDGs ----
test_that("create_query_data Test 6: SDGs", {
  sdg <- query(
    prefix = "SDG01",
    id = auto,
    definition = basket_select(
      name = "5-aminosalicylates for ulcerative colitis",
      scope = NA_character_,
      type = "sdg"
    )
  )

  actual_output <- create_query_data(
    queries = list(sdg),
    version = "2019_09",
    get_terms_fun = get_sdg
  )

  expected_output <-
    get_sdg(
      basket_select(
        name = "5-aminosalicylates for ulcerative colitis",
        scope = NA_character_,
        type = "sdg"
      ),
      version = "2019_09"
    ) %>%
    mutate(
      GRPID = 42,
      PREFIX = "SDG01",
      SCOPE = NA_character_,
      VERSION = "2019_09"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

## Test 7: issues error if SDGs without meddra_version are requested ----
test_that("create_query_data Test 7: issues error if SDGs without meddra_version are requested", {
  sdg <- query(
    prefix = "SDG01",
    definition = basket_select(
      name = "5-aminosalicylates for ulcerative colitis",
      scope = NA_character_,
      type = "sdg"
    )
  )

  expect_snapshot(
    create_query_data(
      queries = list(sdg),
      get_terms_fun = get_sdg
    ),
    error = TRUE
  )
})

# query ----
## Test 8: error if name = auto for non SMQs/SDGs ----
test_that("query Test 8: error if name = auto for non SMQs/SDGs", {
  expect_snapshot(
    sdg <- query(
      prefix = "CQ01",
      definition = cqterms
    ),
    error = TRUE
  )
})

## Test 9: error if id = auto for non SMQs/SDGs ----
test_that("query Test 9: error if id = auto for non SMQs/SDGs", {
  expect_snapshot(
    sdg <- query(
      name = "My CQ",
      id = auto,
      prefix = "CQ01",
      definition = cqterms
    ),
    error = TRUE
  )
})

## Test 10: error if invalid definition ----
test_that("query Test 10: error if invalid definition", {
  expect_snapshot(
    sdg <- query(
      name = "My CQ",
      prefix = "CQ01",
      definition = 1
    ),
    error = TRUE
  )
})

# assert_terms ----
## Test 11: error if SRCVAR missing ----
test_that("assert_terms Test 11: error if SRCVAR missing", {
  expect_snapshot(
    assert_terms(
      terms = select(cqterms, -SRCVAR),
      source_text = "my test data"
    ),
    error = TRUE
  )
})

## Test 12: error if SRCVAR and GRPNAME missing ----
test_that("assert_terms Test 12: error if SRCVAR and GRPNAME missing", {
  expect_snapshot(
    assert_terms(
      terms = select(cqterms, -SRCVAR),
      source_text = "my test data",
      expect_grpname = TRUE
    ),
    error = TRUE
  )
})

## Test 13: error if TERMCHAR and TERMNUM missing ----
test_that("assert_terms Test 13: error if TERMCHAR and TERMNUM missing", {
  expect_snapshot(
    assert_terms(
      terms = select(cqterms, SRCVAR),
      source_text = "my test data"
    ),
    error = TRUE
  )
})

## Test 14: error if no data frame ----
test_that("assert_terms Test 14: error if no data frame", {
  expect_snapshot(
    assert_terms(
      terms = 42,
      source_text = "object returned by calling get_mysmq()"
    ),
    error = TRUE
  )
})

## Test 15: error if no observations ----
test_that("assert_terms Test 15: error if no observations", {
  expect_snapshot(
    assert_terms(
      terms = filter(cqterms, TERMNUM == 42),
      source_text = "object returned by calling get_my_smq"
    ),
    error = TRUE
  )
})

## Test 16: error if GRPNAME is missing ----
test_that("assert_terms Test 16: error if GRPNAME is missing", {
  expect_snapshot(
    assert_terms(
      terms = cqterms,
      expect_grpname = TRUE,
      source_text = "object returned by calling get_my_smq"
    ),
    error = TRUE
  )
})

## Test 17: error if GRPID is missing ----
test_that("assert_terms Test 17: error if GRPID is missing", {
  expect_snapshot(
    assert_terms(
      terms = cqterms,
      expect_grpid = TRUE,
      source_text = "object returned by calling get_my_smq"
    ),
    error = TRUE
  )
})

# basket_select ----
## Test 18: error if name and id specified ----
test_that("basket_select Test 18: error if name and id specified", {
  expect_snapshot(
    basket_select(
      name = "My SMQ",
      id = 42,
      scope = "NARROW",
      type = "smq"
    ),
    error = TRUE
  )
})

## Test 19: error if neither name nor id specified ----
test_that("basket_select Test 19: error if neither name nor id specified", {
  expect_snapshot(
    basket_select(scope = "NARROW", type = "smq"),
    error = TRUE
  )
})

## Test 20: error if type is not specified ----
test_that("basket_select Test 20: error if type is not specified", {
  expect_snapshot(
    basket_select(id = 42, scope = "NARROW"),
    error = TRUE
  )
})

# basket_select customized query defined by SMQs extra arguments ----
get_smq_oth <- function(basket_select,
                        version,
                        keep_id = FALSE,
                        temp_env) {
  if (basket_select$scope == "NARROW") {
    end <- 1
  } else {
    end <- 2
  }

  if (is.null(basket_select$name)) {
    basket_select$name <- paste("SMQ name of", basket_select$id)
  }
  terms <- tibble(TERMCHAR = paste(basket_select$name, "Term", c(1:end)))
  terms <- mutate(terms,
    SRCVAR = "AEDECOD",
    GRPNAME = basket_select$name,
    TEST1_VAR = basket_select$TEST1_VAR,
    TEST2_VAR = basket_select$TEST2_VAR
  )
  if (keep_id) {
    mutate(terms, GRPID = 42)
  } else {
    terms
  }
}


## Test 21: basket_select customized query defined by SMQs extra arguments ----
test_that("basket_select Test 21: basket_select customized query defined by SMQs extra arguments", {
  cq <- query(
    prefix = "CQ02",
    name = "Immune-Mediated Meningoencephalitis",
    definition = list(
      basket_select(
        name = "Noninfectious meningitis",
        scope = "NARROW",
        type = "smq",
        TEST1_VAR = "CHECK 1",
        TEST2_VAR = "CHECK 3"
      ),
      basket_select(
        name = "Noninfectious encephalitis",
        scope = "BROAD",
        type = "smq",
        TEST1_VAR = "CHECK 2",
        TEST2_VAR = "CHECK 4"
      )
    )
  )

  actual_output <- create_query_data(
    queries = list(cq),
    version = "20.0",
    get_terms_fun = get_smq_oth
  )

  expected_output <-
    tribble(
      ~TERMCHAR,                           ~TEST1_VAR, ~TEST2_VAR,
      "Noninfectious meningitis Term 1",    "CHECK 1",  "CHECK 3",
      "Noninfectious encephalitis Term 1",  "CHECK 2",  "CHECK 4",
      "Noninfectious encephalitis Term 2",  "CHECK 2",  "CHECK 4",
    ) %>%
    mutate(
      SRCVAR = "AEDECOD",
      GRPNAME = "Immune-Mediated Meningoencephalitis",
      PREFIX = "CQ02",
      VERSION = "20.0"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("PREFIX", "TERMCHAR")
  )
})

# format.basket_select ----
## Test 22: formatting is correct (id specified) ----
test_that("format.basket_select Test 22: formatting is correct (id specified)", {
  expect_equal(
    format(basket_select(
      id = 42,
      scope = "NARROW",
      type = "smq",
      newvar = 1
    )),
    "basket_select(name = NULL, id = 42, scope = \"NARROW\", type = \"smq\", newvar = 1)"
  )
})

## Test 23: formatting is correct (name specified) ----
test_that("format.basket_select Test 23: formatting is correct (name specified)", {
  expect_equal(
    format(basket_select(name = "My SDG", type = "sdg", scope = NA_character_)),
    "basket_select(name = \"My SDG\", id = NULL, scope = \"NA\", type = \"sdg\")"
  )
})
