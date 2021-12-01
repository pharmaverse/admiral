context("test-read_query_data")

get_smq <- function(query,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  if (query$scope == "NARROW") {
    end <- 3
  }
  else {
    end <- 5
  }

  terms <- tibble(TERM_NAME = paste(query$name, "Term", c(1:end), "(", version, ")"))
  terms <- mutate(terms, TERM_LEVEL = "AEDECOD", QUERY_NAME = query$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  }
  else {
    terms
  }
}

get_sdg <- function(query,
                    version,
                    keep_id = FALSE,
                    temp_env) {
  terms <- tibble(TERM_NAME = paste(query$name, "Term", c(1:end)))
  terms <- mutate(terms, TERM_LEVEL = "CMDECOD", QUERY_NAME = query$name)
  if (keep_id) {
    mutate(terms, QUERY_ID = 42)
  }
  else {
    terms
  }
}

cqterms <- tribble(
  ~TERM_NAME, ~TERM_ID,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L) %>%
  mutate(TERM_LEVEL = "AEDECOD")

test_that("customized query defined by terms", {
  cq <- query(prefix = "CQ01",
              name = "Application Site Issues",
              definition = cqterms)

  actual_output <- read_query_data(queries = list(cq))

  expected_output <- cqterms %>% mutate(QUERY_NAME = "Application Site Issues",
                                        VAR_PREFIX = "CQ01")

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})

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

  actual_output <- read_query_data(queries = list(cq),
                                   meddra_version = "20.0",
                                   get_smq_fun = get_smq)

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

  actual_output <- read_query_data(queries = list(cq),
                                   meddra_version = "20.1",
                                   get_smq_fun = get_smq)

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

test_that("SMQs", {
  pregsqm <- query(
    prefix = "SMQ02",
    id = auto,
    scope = "NARROW",
    scope_num = 2,
    definition = smq_select(name = "Pregnancy and neonatal topics (SMQ)")
  )

  pneuaegt <- query(prefix = "SMQ04",
                    scope = "BROAD",
                    definition = smq_select(id = 8050L))

  actual_output <- read_query_data(queries = list(pregsqm, pneuaegt),
                                   meddra_version = "20.0",
                                   get_smq_fun = get_smq)

  expected_output <-
    bind_rows(
      get_smq(
        smq_select(name = "Pregnancy and neonatal topics (SMQ)",
                   scope =  "NARROW"),
        version = "20.0"
      ) %>%
        mutate(QUERY_NAME = "Pregnancy and neonatal topics (SMQ)",
               VAR_PREFIX = "SMQ02")            ,
      get_smq(
        smq_select(id = 8050L,
                   scope = "BROAD"),
        version = "20.0"
      ) %>%
        mutate(VAR_PREFIX = "SMQ04")
    )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("VAR_PREFIX", "TERM_NAME")
  )
})
