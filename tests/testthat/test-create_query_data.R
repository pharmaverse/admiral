
get_smq <- function(smq_select,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  if (smq_select$scope == "NARROW") {
    end <- 3
  }
  else {
    end <- 5
  }

  if (is.null(smq_select$name)) {
    smq_select$name <- paste("SMQ name of", smq_select$id)
  }
  terms <- tibble(TERM_NAME = paste(smq_select$name, "Term", c(1:end), "(", version, ")"))
  terms <- mutate(terms, TERM_LEVEL = "AEDECOD", QUERY_NAME = smq_select$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  }
  else {
    terms
  }
}

get_sdg <- function(sdg_select,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  terms <- tibble(TERM_NAME = paste(sdg_select$name, "Term", c(1:4)))
  terms <- mutate(terms, TERM_LEVEL = "CMDECOD", QUERY_NAME = sdg_select$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  }
  else {
    terms
  }
}

cqterms <- tibble::tribble(
  ~TERM_NAME, ~TERM_ID,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L) %>%
  mutate(TERM_LEVEL = "AEDECOD")

# customized query defined by terms ----
test_that("customized query defined by terms", {
  cq <- query(prefix = "CQ01",
              name = "Application Site Issues",
              definition = cqterms)

  actual_output <- create_query_data(queries = list(cq))

  expected_output <- cqterms %>% mutate(QUERY_NAME = "Application Site Issues",
                                        VAR_PREFIX = "CQ01")

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# customized query defined by SMQs ----
test_that("customized query defined by SMQs", {
  cq <- query(
    prefix = "CQ02",
    name = "Immune-Mediated Meningoencephalitis",
    definition = list(
      smq_select(
        name = "Noninfectious meningitis",
        scope =  "NARROW"
      ),
      smq_select(
        name = "Noninfectious encephalitis",
        scope = "BROAD"
      )))

  actual_output <- create_query_data(
    queries = list(cq),
    meddra_version = "20.0",
    get_smq_fun = get_smq
  )

  expected_output <-
    bind_rows(get_smq(smq_select(name = "Noninfectious meningitis",
                                 scope =  "NARROW"),
                      version = "20.0"),
              get_smq(smq_select(name = "Noninfectious encephalitis",
                                 scope = "BROAD"),
                      version = "20.0")) %>%
    mutate(QUERY_NAME = "Immune-Mediated Meningoencephalitis",
           VAR_PREFIX = "CQ02")

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

test_that("customized query defined by terms and SMQs", {
  cq <- query(
    prefix = "CQ03",
    name = "Immune-Mediated Meningoencephalitis or Application Site Issues",
    definition = list(
      smq_select(
        name = "Noninfectious meningitis",
        scope =  "NARROW"
      ),
      cqterms,
      smq_select(
        name = "Noninfectious encephalitis",
        scope = "BROAD"
      )))

  actual_output <- create_query_data(
    queries = list(cq),
    meddra_version = "20.1",
    get_smq_fun = get_smq
  )

  expected_output <-
    bind_rows(get_smq(smq_select(name = "Noninfectious meningitis",
                                 scope =  "NARROW"),
                      version = "20.1"),
              cqterms,
              get_smq(smq_select(name = "Noninfectious encephalitis",
                                 scope = "BROAD"),
                      version = "20.1")) %>%
    mutate(QUERY_NAME = "Immune-Mediated Meningoencephalitis or Application Site Issues",
           VAR_PREFIX = "CQ03")

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# SMQs ----
test_that("SMQs", {
  pregsmq <- query(
    prefix = "SMQ02",
    id = 13,
    add_scope_num = TRUE,
    definition = smq_select(name = "Pregnancy and neonatal topics (SMQ)",
                            scope = "NARROW")
  )

  pneuaegt <- query(prefix = "SMQ04",
                    definition = smq_select(id = 8050L,
                                            scope = "BROAD"))

  actual_output <-
    create_query_data(
      queries = list(pregsmq, pneuaegt),
      meddra_version = "20.0",
      get_smq_fun = get_smq
    )

  expected_output <-
    bind_rows(
      get_smq(
        smq_select(name = "Pregnancy and neonatal topics (SMQ)",
                   scope = "NARROW"),
        version = "20.0"
      ) %>%
        mutate(QUERY_NAME = "Pregnancy and neonatal topics (SMQ)",
               QUERY_ID = 13,
               QUERY_SCOPE = "NARROW",
               QUERY_SCOPE_NUM = 2,
               VAR_PREFIX = "SMQ02"),
      get_smq(smq_select(id = 8050L,
                         scope = "BROAD"),
              version = "20.0") %>%
        mutate(QUERY_SCOPE = "BROAD",
               VAR_PREFIX = "SMQ04")
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

# issues error if SMQs without get_smq_fun are requested ----
test_that("issues error if SMQs without get_smq_fun are requested", {
  pregsmq <- query(
    prefix = "SMQ02",
    definition = smq_select(name = "Pregnancy and neonatal topics (SMQ)",
                            scope = "NARROW")
  )

  expect_error(
    create_query_data(queries = list(pregsmq),
                      meddra_version = "20.0"),
    regexp = "^get_smq_fun is not specified. This is expected for SMQs.*")
})

# issues error if SMQs without meddra_version are requested ----
test_that("issues error if SMQs without meddra_version are requested", {
  pregsmq <- query(
    prefix = "SMQ02",
    definition = smq_select(name = "Pregnancy and neonatal topics (SMQ)",
                            scope = "NARROW")
  )

  expect_error(
    create_query_data(queries = list(pregsmq),
                      get_smq_fun = get_smq),
    regexp = "^meddra_version is not specified. This is expected for SMQs.*")
})

# SDGs ----
test_that("SDGs", {
  sdg <- query(
    prefix = "SDG01",
    id = auto,
    definition = sdg_select(name = "5-aminosalicylates for ulcerative colitis")
  )

  actual_output <- create_query_data(queries = list(sdg),
                                   whodd_version = "2019_09",
                                   get_sdg_fun = get_sdg)

  expected_output <-
    get_sdg(sdg_select(name = "5-aminosalicylates for ulcerative colitis"),
            version = "2019_09") %>%
    mutate(
      QUERY_ID = 42,
      VAR_PREFIX = "SDG01")

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})
# issues error if SDGs without get_sdg_fun are requested ----
test_that("issues error if SDGs without get_sdg_fun are requested", {
  sdg <- query(
    prefix = "SDG01",
    definition = sdg_select(name = "5-aminosalicylates for ulcerative colitis")
  )

  expect_error(
    create_query_data(queries = list(sdg),
                      whodd_version = "2019_09"),
    regexp = "^get_sdg_fun is not specified. This is expected for SDGs.*")
})

# issues error if SDGs without meddra_version are requested ----
test_that("issues error if SDGs without meddra_version are requested", {
  sdg <- query(
    prefix = "SDG01",
    definition = sdg_select(name = "5-aminosalicylates for ulcerative colitis")
  )

  expect_error(
    create_query_data(queries = list(sdg),
                      get_sdg_fun = get_sdg),
    regexp = "^whodd_version is not specified. This is expected for SDGs.*")
})

# query: error: add_scope_num = TRUE for non SMQs ----
test_that("query: error: add_scope_num = TRUE for non SMQs", {
  expect_error(
    sdg <- query(
      prefix = "SDG01",
      add_scope_num = TRUE,
      definition = sdg_select(name = "5-aminosalicylates for ulcerative colitis")
    ),
    regexp = "`add_scope_num == TRUE` must be used for SMQs only.",
    fixed = TRUE)
})

# query: error: name = auto for non SMQs/SDGs ----
test_that("query: error: name = auto for non SMQs/SDGs", {
  expect_error(
    sdg <- query(
      prefix = "CQ01",
      definition = cqterms
    ),
    regexp = "^The auto keyword can be used for SMQs and SDGs only.*")
})

# query: error: name = id for non SMQs/SDGs ----
test_that("query: error: name = id for non SMQs/SDGs", {
  expect_error(
    sdg <- query(
      name = "My CQ",
      id = auto,
      prefix = "CQ01",
      definition = cqterms
    ),
    regexp = "^The auto keyword can be used for SMQs and SDGs only.*")
})

# query: error: definition is list with non dataframe or smq_select elements ----
test_that("query: error: definition is list with non dataframe or smq_select elements", {
  expect_error(
    sdg <- query(
      name = "My CQ",
      prefix = "CQ01",
      definition = list(sdg_select(name = "5-aminosalicylates for ulcerative colitis"))
    ),
    regexp =
      paste0(
        "^Each element of the list in the definition field must be a data frame or",
        " an object of class `smq_select` but the following are not:.*"
      )
  )
})

# query: error: invalid definition ----
test_that("query: error: invalid definition", {
  expect_error(
    sdg <- query(
      name = "My CQ",
      prefix = "CQ01",
      definition = 1
    ),
    regexp =
      paste0(
        "^`definition` expects a `smq_select` or `sdg_select` object,",
        " a data frame, or a list of data frames and `smq_select` objects.*"
      )
  )
})

# assert_terms: error: TERM_LEVEL missing ----
test_that("assert_terms: error: TERM_LEVEL missing", {
  expect_error(
    assert_terms(terms = select(cqterms, -TERM_LEVEL),
                 source_text = "my test data"),
    regexp = "Required variable `TERM_LEVEL` is missing in my test data.",
    fixed = TRUE)
})

# assert_terms: error: TERM_NAME and TERM_ID missing ----
test_that("assert_terms: error: TERM_NAME and TERM_ID missing", {
  expect_error(
    assert_terms(terms = select(cqterms, TERM_LEVEL),
                 source_text = "my test data"),
    regexp = paste0(
      "Variable `TERM_NAME` or `TERM_ID` is required.\n",
      "None of them is in my test data.\nProvided variables: `TERM_LEVEL`"
    ),
    fixed = TRUE)
})

# assert_terms: error: no data frame ----
test_that("assert_terms: error: no data frame", {
  expect_error(
    assert_terms(terms = 42,
                 source_text = "object returned by calling get_mysmq()"),
    regexp = "object returned by calling get_mysmq() is not a data frame but `42`.",
    fixed = TRUE)
})

# assert_terms: error: no observations ----
test_that("assert_terms: error: no observations", {
  expect_error(
    assert_terms(terms = filter(cqterms, TERM_ID == 42),
                 source_text = "object returned by calling get_my_smq"),
    regexp = "object returned by calling get_my_smq does not contain any observations.",
    fixed = TRUE)
})

# assert_terms: error: QUERY_NAME is missing ----
test_that("assert_terms: error: QUERY_NAME is missing", {
  expect_error(
    assert_terms(terms = cqterms,
                 expect_query_name = TRUE,
                 source_text = "object returned by calling get_my_smq"),
    regexp = "Required variable `QUERY_NAME` is missing in object returned by calling get_my_smq.",
    fixed = TRUE)
})

# assert_terms: error: QUERY_ID is missing ----
test_that("assert_terms: error: QUERY_ID is missing", {
  expect_error(
    assert_terms(terms = cqterms,
                 expect_query_id = TRUE,
                 source_text = "object returned by calling get_my_smq"),
    regexp = "Required variable `QUERY_ID` is missing in object returned by calling get_my_smq.",
    fixed = TRUE)
})

# smq_select: error: name and id specified ----
test_that("smq_select: error: name and id specified", {
  expect_error(
    smq_select(name = "My SMQ",
               id = 42,
               scope = "NARROW"),
    regexp = "Either id or name has to be null.",
    fixed = TRUE)
})

# smq_select: error: neither name nor id specified ----
test_that("smq_select: error: neither name nor id specified", {
  expect_error(
    smq_select(scope = "NARROW"),
    regexp = "Either id or name has to be non null.",
    fixed = TRUE)
})

# sdg_select: error: name and id specified ----
test_that("sdg_select: error: name and id specified", {
  expect_error(
    sdg_select(name = "My SDG",
               id = 42),
    regexp = "Either id or name has to be null.",
    fixed = TRUE)
})

# format.smq_select: formatting is correct ----
test_that("format.smq_select: formatting is correct", {
  expect_equal(format(smq_select(id = 42,
                                 scope = "NARROW")),
               "smq_select(name = NULL, id = 42, scope = \"NARROW\")")
})

# sdg_select: error: neither name nor id specified ----
test_that("sdg_select: error: neither name nor id specified", {
  expect_error(
    sdg_select(),
    regexp = "Either id or name has to be non null.",
    fixed = TRUE)
})

# format.sdg_select: formatting is correct ----
test_that("format.sdg_select: formatting is correct", {
  expect_equal(format(sdg_select(name = "My SDG")),
               "sdg_select(name = \"My SDG\", id = NULL)")
})
