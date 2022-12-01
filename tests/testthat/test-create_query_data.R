
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
  terms <- tibble(TERM_NAME = paste(basket_select$name, "Term", c(1:end), "(", version, ")"))
  terms <- mutate(terms, TERM_LEVEL = "AEDECOD", QUERY_NAME = basket_select$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  } else {
    terms
  }
}

get_sdg <- function(basket_select,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  terms <- tibble(TERM_NAME = paste(basket_select$name, "Term", c(1:4)))
  terms <- mutate(terms, TERM_LEVEL = "CMDECOD", QUERY_NAME = basket_select$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  } else {
    terms
  }
}

cqterms <- tibble::tribble(
  ~TERM_NAME, ~TERM_ID,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L
) %>%
  mutate(TERM_LEVEL = "AEDECOD")
# create_query_data ----
# customized query defined by terms ----
## Test 1: customized query defined by terms ----
test_that("create_query_data Test 1: customized query defined by terms", {
  cq <- query(
    prefix = "CQ01",
    name = "Application Site Issues",
    definition = cqterms
  )

  actual_output <- create_query_data(queries = list(cq))

  expected_output <- cqterms %>% mutate(
    QUERY_NAME = "Application Site Issues",
    VAR_PREFIX = "CQ01"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# customized query defined by SMQs ----
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
      QUERY_NAME = "Immune-Mediated Meningoencephalitis",
      VAR_PREFIX = "CQ02",
      VERSION = "20.0"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
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
      QUERY_NAME = "Immune-Mediated Meningoencephalitis or Application Site Issues",
      VAR_PREFIX = "CQ03",
      VERSION = "20.1"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# SMQs ----
## Test 4: SMQs ----
test_that("SMQs Test 4: SMQs", {
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
          QUERY_NAME = "Pregnancy and neonatal topics (SMQ)",
          QUERY_ID = 13,
          QUERY_SCOPE = "NARROW",
          QUERY_SCOPE_NUM = 2,
          VAR_PREFIX = "SMQ02"
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
          QUERY_SCOPE = "BROAD",
          VAR_PREFIX = "SMQ04"
        )
    ) %>%
    mutate(
      VERSION = "20.0"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# issues error if SMQs without meddra_version are requested ----
## Test 5: issues error if SMQs without meddra_version are requested ----
test_that("SMQs Test 5: issues error if SMQs without meddra_version are requested", {
  pregsmq <- query(
    prefix = "SMQ02",
    definition = basket_select(
      name = "Pregnancy and neonatal topics (SMQ)",
      scope = "NARROW",
      type = "smq"
    )
  )

  expect_error(
    create_query_data(
      queries = list(pregsmq),
      get_terms_fun = get_smq
    ),
    regexp = "^version is not specified. This is expected for baskets.*"
  )
})

# SDGs ----
## Test 6: SDGs ----
test_that("SDGs Test 6: SDGs", {
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
      QUERY_ID = 42,
      VAR_PREFIX = "SDG01",
      QUERY_SCOPE = NA_character_,
      VERSION = "2019_09"
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# issues error if SDGs without meddra_version are requested ----
## Test 7: issues error if SDGs without meddra_version are requested ----
test_that("SDGs Test 7: issues error if SDGs without meddra_version are requested", {
  sdg <- query(
    prefix = "SDG01",
    definition = basket_select(
      name = "5-aminosalicylates for ulcerative colitis",
      scope = NA_character_,
      type = "sdg"
    )
  )

  expect_error(
    create_query_data(
      queries = list(sdg),
      get_terms_fun = get_sdg
    ),
    regexp = "^version is not specified. This is expected for baskets.*"
  )
})

# query() ----
# query: error: name = auto for non SMQs/SDGs ----
## Test 8: query: error: name = auto for non SMQs/SDGs ----
test_that("SDGs Test 8: query: error: name = auto for non SMQs/SDGs", {
  expect_error(
    sdg <- query(
      prefix = "CQ01",
      definition = cqterms
    ),
    regexp = "^The auto keyword can be used for baskets only.*"
  )
})

# query: error: name = id for non SMQs/SDGs ----
## Test 9: query: error: name = id for non SMQs/SDGs ----
test_that("SDGs Test 9: query: error: name = id for non SMQs/SDGs", {
  expect_error(
    sdg <- query(
      name = "My CQ",
      id = auto,
      prefix = "CQ01",
      definition = cqterms
    ),
    regexp = "^The auto keyword can be used for baskets only.*"
  )
})

# query: error: invalid definition ----
## Test 10: query: error: invalid definition ----
test_that("SDGs Test 10: query: error: invalid definition", {
  expect_error(
    sdg <- query(
      name = "My CQ",
      prefix = "CQ01",
      definition = 1
    ),
    regexp =
      paste0(
        "^`definition` expects a `basket_select` object,",
        " a data frame, or a list of data frames and `basket_select` objects*"
      )
  )
})

# assert_terms ----
# assert_terms: error: TERM_LEVEL missing ----
## Test 11: assert_terms: error: TERM_LEVEL missing ----
test_that("assert_terms Test 11: assert_terms: error: TERM_LEVEL missing", {
  expect_error(
    assert_terms(
      terms = select(cqterms, -TERM_LEVEL),
      source_text = "my test data"
    ),
    regexp = "Required variable `TERM_LEVEL` is missing in my test data.",
    fixed = TRUE
  )
})

# assert_terms: error: TERM_NAME and TERM_ID missing ----
## Test 12: assert_terms: error: TERM_NAME and TERM_ID missing ----
test_that("assert_terms Test 12: assert_terms: error: TERM_NAME and TERM_ID missing", {
  expect_error(
    assert_terms(
      terms = select(cqterms, TERM_LEVEL),
      source_text = "my test data"
    ),
    regexp = paste0(
      "Variable `TERM_NAME` or `TERM_ID` is required.\n",
      "None of them is in my test data.\nProvided variables: `TERM_LEVEL`"
    ),
    fixed = TRUE
  )
})

# assert_terms: error: no data frame ----
## Test 13: assert_terms: error: no data frame ----
test_that("assert_terms Test 13: assert_terms: error: no data frame", {
  expect_error(
    assert_terms(
      terms = 42,
      source_text = "object returned by calling get_mysmq()"
    ),
    regexp = "object returned by calling get_mysmq() is not a data frame but `42`.",
    fixed = TRUE
  )
})

# assert_terms: error: no observations ----
## Test 14: assert_terms: error: no observations ----
test_that("assert_terms Test 14: assert_terms: error: no observations", {
  expect_error(
    assert_terms(
      terms = filter(cqterms, TERM_ID == 42),
      source_text = "object returned by calling get_my_smq"
    ),
    regexp = "object returned by calling get_my_smq does not contain any observations.",
    fixed = TRUE
  )
})

# assert_terms: error: QUERY_NAME is missing ----
## Test 15: assert_terms: error: QUERY_NAME is missing ----
test_that("assert_terms Test 15: assert_terms: error: QUERY_NAME is missing", {
  expect_error(
    assert_terms(
      terms = cqterms,
      expect_query_name = TRUE,
      source_text = "object returned by calling get_my_smq"
    ),
    regexp = "Required variable `QUERY_NAME` is missing in object returned by calling get_my_smq.",
    fixed = TRUE
  )
})

# assert_terms: error: QUERY_ID is missing ----
## Test 16: assert_terms: error: QUERY_ID is missing ----
test_that("assert_terms Test 16: assert_terms: error: QUERY_ID is missing", {
  expect_error(
    assert_terms(
      terms = cqterms,
      expect_query_id = TRUE,
      source_text = "object returned by calling get_my_smq"
    ),
    regexp = "Required variable `QUERY_ID` is missing in object returned by calling get_my_smq.",
    fixed = TRUE
  )
})

# basket_select ----
# basket_select: error: name and id specified ----
## Test 17: basket_select: error: name and id specified ----
test_that("basket_select Test 17: basket_select: error: name and id specified", {
  expect_error(
    basket_select(
      name = "My SMQ",
      id = 42,
      scope = "NARROW",
      type = "smq"
    ),
    regexp = "Either id or name has to be null.",
    fixed = TRUE
  )
})

# basket_select: error: neither name nor id specified ----
## Test 18: basket_select: error: neither name nor id specified ----
test_that("basket_select Test 18: basket_select: error: neither name nor id specified", {
  expect_error(
    basket_select(scope = "NARROW", type = "smq"),
    regexp = "Either id or name has to be non null.",
    fixed = TRUE
  )
})

# basket_select: error: name and id specified ----
## Test 19: basket_select: error: name and id specified ----
test_that("basket_select Test 19: basket_select: error: name and id specified", {
  expect_error(
    basket_select(
      name = "My SDG",
      id = 42,
      scope = NA_character_,
      type = "sdg"
    ),
    regexp = "Either id or name has to be null.",
    fixed = TRUE
  )
})

# format.basket_select: formatting is correct ----
## Test 20: format.basket_select: formatting is correct ----
test_that("basket_select Test 20: format.basket_select: formatting is correct", {
  expect_equal(
    format(basket_select(
      id = 42,
      scope = "NARROW",
      type = "smq"
    )),
    "basket_select(name = NULL, id = 42, scope = \"NARROW\", type = \"smq\")"
  )
})

# basket_select: error: neither name nor id specified ----
## Test 21: basket_select: error: neither name nor id specified ----
test_that("basket_select Test 21: basket_select: error: neither name nor id specified", {
  expect_error(
    basket_select(type = "sdg", scope = NA_character_),
    regexp = "Either id or name has to be non null.",
    fixed = TRUE
  )
})

# basket_select: error: type is not specified ----
## Test 22: basket_select: error: type is not specified ----
test_that("basket_select Test 22: basket_select: error: type is not specified", {
  expect_error(
    basket_select(id = 42, scope = "NARROW"),
    regexp = "argument \"type\" is missing, with no default",
    fixed = TRUE
  )
})

# format.basket_select ----
# format.basket_select: formatting is correct ----
## Test 23: format.basket_select: formatting is correct ----
test_that("format.basket_select Test 23: format.basket_select: formatting is correct", {
  expect_equal(
    format(basket_select(name = "My SDG", type = "sdg", scope = NA_character_)),
    "basket_select(name = \"My SDG\", id = NULL, scope = \"NA\", type = \"sdg\")"
  )
})
