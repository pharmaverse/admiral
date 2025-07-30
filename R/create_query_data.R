#' Creates a queries dataset as input dataset to the `dataset_queries` argument in
#' `derive_vars_query()`
#'
#' Creates a queries dataset as input dataset to the `dataset_queries` argument
#' in the `derive_vars_query()` function as defined in the
#' `vignette("queries_dataset")`.
#'
#' @param queries List of queries
#'
#'   A list of `query()` objects is expected.
#' @param version Dictionary version
#'
#'   The dictionary version used for coding the terms should be specified.
#'   If any of the queries is a basket (SMQ, SDG, ....) or a customized query
#'   including a basket, the parameter needs to be specified.
#'
#' @permitted A character string (the expected format is
#'   company-specific)
#'
#' @param get_terms_fun Function which returns the terms
#'
#'   For each query specified for the `queries` parameter referring to a basket
#'   (i.e., those where the `definition` field is set to a `basket_select()`
#'   object or a list which contains at least one `basket_select()` object) the
#'   specified function is called to retrieve the terms defining the query.
#'   This function is not provided by admiral as it is company specific, i.e.,
#'   it has to be implemented at company level.
#'
#'   The function must return a dataset with all the terms defining the basket.
#'   The output dataset must contain the following variables.
#'
#'   - `SRCVAR`: the variable to be used for defining a term of the basket,
#'    e.g., `AEDECOD`
#'   - `TERMCHAR`: the name of the term if the variable `SRCVAR` is
#'   referring to is character
#'   - `TERMNUM` the numeric id of the term if the variable `SRCVAR` is
#'   referring to is numeric
#'   - `GRPNAME`: the name of the basket. The values must be the same for
#'   all observations.
#'
#'   The function must provide the following parameters
#'
#'   - `basket_select`: A `basket_select()` object.
#'   - `version`: The dictionary version. The value specified for the
#'   `version` in the `create_query_data()` call is passed to this
#'   parameter.
#'   - `keep_id`: If set to `TRUE`, the output dataset must contain the
#'   `GRPID` variable. The variable must be set to the numeric id of the basket.
#'   - `temp_env`: A temporary environment is passed to this parameter. It can
#'   be used to store data which is used for all baskets in the
#'   `create_query_data()` call. For example if SMQs need to be read from a
#'   database all SMQs can be read and stored in the environment when the first
#'   SMQ is handled. For the other SMQs the terms can be retrieved from the
#'   environment instead of accessing the database again.
#'
#' @details
#'
#'   For each `query()` object listed in the `queries` argument, the terms belonging
#'   to the query (`SRCVAR`, `TERMCHAR`, `TERMNUM`) are determined with respect
#'   to the `definition` field of the query: if the definition field of the
#'   `query()` object is
#'
#'   * a `basket_select()` object, the terms are read from the basket
#'   database by calling the function specified for the `get_terms_fun` parameter.
#'   * a data frame, the terms stored in the data frame are used.
#'   * a list of data frames and `basket_select()` objects, all terms from
#'   the data frames and all terms read from the basket database referenced by the
#'   `basket_select()` objects are collated.
#'
#' The following variables (as described in `vignette("queries_dataset")` are
#' created:
#'
#'   * `PREFIX`: Prefix of the variables to be created by
#'   `derive_vars_query()` as specified by the `prefix` element.
#'   * `GRPNAME`: Name of the query as specified by the `name` element.
#'   * `GRPID`: Id of the query as specified by the `id` element. If the `id`
#'   element is not specified for a query, the variable is set to `NA`. If the
#'   `id` element is not specified for any query, the variable is not created.
#'   * `SCOPE`: scope of the query as specified by the `scope` element of
#'   the `basket_select()` object. For queries not defined by a `basket_select()`
#'   object, the variable is set to `NA`. If none of the queries is defined by a
#'   `basket_select()` object, the variable is not created.
#'   * `SCOPEN`: numeric scope of the query. It is set to `1` if the
#'   scope is broad. Otherwise it is set to `2`. If the `add_scope_num` element
#'   equals `FALSE`, the variable is set to `NA`. If the `add_scope_num` element
#'   equals `FALSE` for all baskets or none of the queries is an basket , the variable
#'   is not created.
#'   * `SRCVAR`: Name of the variable used to identify the terms.
#'   * `TERMCHAR`: Value of the term variable if it is a character variable.
#'   * `TERMNUM`: Value of the term variable if it is a numeric variable.
#'   * `VERSION`: Set to the value of the `version` argument. If it is not
#'   specified, the variable is not created.
#'
#'
#' @return A dataset to be used as input dataset to the `dataset_queries`
#'   argument in `derive_vars_query()`
#'
#' @family create_aux
#' @keywords create_aux
#'
#' @seealso [derive_vars_query()], [query()], [basket_select()],
#'   `vignette("queries_dataset")`
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(pharmaversesdtm)
#' library(admiral)
#'
#' # creating a query dataset for a customized query
#' cqterms <- tribble(
#'   ~TERMCHAR, ~TERMNUM,
#'   "APPLICATION SITE ERYTHEMA", 10003041L,
#'   "APPLICATION SITE PRURITUS", 10003053L
#' ) %>%
#'   mutate(SRCVAR = "AEDECOD")
#'
#' cq <- query(
#'   prefix = "CQ01",
#'   name = "Application Site Issues",
#'   definition = cqterms
#' )
#'
#' create_query_data(queries = list(cq))
#'
#' # create a query dataset for SMQs
#' pregsmq <- query(
#'   prefix = "SMQ02",
#'   id = auto,
#'   definition = basket_select(
#'     name = "Pregnancy and neonatal topics (SMQ)",
#'     scope = "NARROW",
#'     type = "smq"
#'   )
#' )
#'
#' bilismq <- query(
#'   prefix = "SMQ04",
#'   definition = basket_select(
#'     id = 20000121L,
#'     scope = "BROAD",
#'     type = "smq"
#'   )
#' )
#'
#' # The get_terms function from pharmaversesdtm is used for this example.
#' # In a real application a company-specific function must be used.
#' create_query_data(
#'   queries = list(pregsmq, bilismq),
#'   get_terms_fun = pharmaversesdtm:::get_terms,
#'   version = "20.1"
#' )
#'
#' # create a query dataset for SDGs
#' sdg <- query(
#'   prefix = "SDG01",
#'   id = auto,
#'   definition = basket_select(
#'     name = "5-aminosalicylates for ulcerative colitis",
#'     scope = NA_character_,
#'     type = "sdg"
#'   )
#' )
#'
#' # The get_terms function from pharmaversesdtm is used for this example.
#' # In a real application a company-specific function must be used.
#' create_query_data(
#'   queries = list(sdg),
#'   get_terms_fun = pharmaversesdtm:::get_terms,
#'   version = "2019-09"
#' )
#'
#' # creating a query dataset for a customized query including SMQs
#' # The get_terms function from pharmaversesdtm is used for this example.
#' # In a real application a company-specific function must be used.
#' create_query_data(
#'   queries = list(
#'     query(
#'       prefix = "CQ03",
#'       name = "Special issues of interest",
#'       definition = list(
#'         basket_select(
#'           name = "Pregnancy and neonatal topics (SMQ)",
#'           scope = "NARROW",
#'           type = "smq"
#'         ),
#'         cqterms
#'       )
#'     )
#'   ),
#'   get_terms_fun = pharmaversesdtm:::get_terms,
#'   version = "20.1"
#' )
create_query_data <- function(queries,
                              version = NULL,
                              get_terms_fun = NULL) {
  # check parameters
  assert_character_scalar(version, optional = TRUE)
  assert_function(get_terms_fun,
    params = c("basket_select", "version", "keep_id", "temp_env"),
    optional = TRUE
  )

  walk(queries, validate_query)

  # read queries
  temp_env <- new.env(parent = emptyenv())
  query_data <- vector("list", length(queries))
  for (i in seq_along(queries)) {
    # get term names and term variable
    if (inherits(queries[[i]]$definition, "basket_select")) {
      # query is a basket
      query_data[[i]] <- get_terms_from_db(
        version = version,
        fun = get_terms_fun,
        fun_name = deparse(substitute(get_terms_fun)),
        queries = queries,
        definition = queries[[i]]$definition,
        expect_grpname = TRUE,
        expect_grpid = !is.null(queries[[i]]$id),
        i = i,
        temp_env = temp_env
      )
      query_data[[i]] <- mutate(query_data[[i]],
        SCOPE = queries[[i]]$definition$scope
      )
      if (queries[[i]]$add_scope_num) {
        query_data[[i]] <-
          mutate(query_data[[i]],
            SCOPEN = if_else(SCOPE == "BROAD", 1, 2)
          )
      }
    } else if (is.data.frame(queries[[i]]$definition)) {
      # query is a customized query
      query_data[[i]] <- queries[[i]]$definition
    } else if (is.list(queries[[i]]$definition)) {
      # query is defined by customized queries and baskets
      definition <- queries[[i]]$definition
      terms <- vector("list", length(definition))
      for (j in seq_along(definition)) {
        if (is.data.frame(definition[[j]])) {
          terms[[j]] <- definition[[j]]
        } else {
          terms[[j]] <- get_terms_from_db(
            version = version,
            fun = get_terms_fun,
            fun_name = deparse(substitute(get_terms_fun)),
            queries = queries,
            definition = definition[[j]],
            i = i,
            temp_env = temp_env
          )
        }
      }
      query_data[[i]] <- bind_rows(terms)
    }

    # add mandatory variables
    query_data[[i]] <- mutate(
      query_data[[i]],
      PREFIX = queries[[i]]$prefix
    )

    if (!is_auto(queries[[i]]$name)) {
      query_data[[i]] <- mutate(
        query_data[[i]],
        GRPNAME = queries[[i]]$name
      )
    }

    # add optional variables
    if (!is.null(queries[[i]]$id) && !is_auto(queries[[i]]$id)) {
      query_data[[i]] <- mutate(query_data[[i]],
        GRPID = queries[[i]]$id
      )
    }
  }
  queries <- bind_rows(query_data)

  if (!is.null(version)) {
    queries <- mutate(
      queries,
      VERSION = version
    )
  }
  queries
}

#' Get Terms from the Queries Database
#'
#' The function checks if all requirements to access the database are fulfilled
#' (version and access function are available, see `assert_db_requirements()`),
#' reads the terms from the database, and checks if the dataset with the terms
#' is in the expected format (see `assert_terms()`).
#'
#' @param version Version
#'
#'   The version must be non null. Otherwise, an error is issued. The value is
#'   passed to the access function (`fun`).
#'
#' @param fun Access function
#'
#'   The access function must be non null. Otherwise, an error is issued. The
#'   function is called to retrieve the terms.
#'
#' @param fun_name Name of access function
#'
#'   The character name of the access function, usually created with
#'   `deparse(substitute(fun))`. This must be non null. Otherwise, an error is issued.
#'
#' @param queries Queries
#'
#'   List of all queries passed to `create_query_data()`. It is used for error
#'   messages.
#'
#' @param definition Definition of the query
#'
#'   The definition is passed to the access function. It defines which terms are
#'   returned.
#'
#' @param expect_grpname Is `GRPNAME` expected in the output dataset?
#'
#' @param expect_grpid Is `GRPID` expected in the output dataset?
#'
#' @param i Index of `definition` in `queries`
#'
#'   The value is used for error messages.
#'
#' @param temp_env Temporary environment
#'
#'   The value is passed to the access function.
#'
#' @keywords internal
#'
#' @return Output dataset of the access function
#'
get_terms_from_db <- function(version,
                              fun,
                              fun_name,
                              queries,
                              definition,
                              expect_grpname = FALSE,
                              expect_grpid = FALSE,
                              i,
                              temp_env) {
  assert_character_scalar(fun_name)

  assert_db_requirements(
    version = version,
    version_arg_name = deparse(substitute(version)),
    fun = fun,
    fun_arg_name = deparse(substitute(fun)),
    queries = queries,
    i = i
  )

  terms <- tryCatch(
    fun(
      version = version,
      basket_select = definition,
      keep_id = expect_grpid,
      temp_env = temp_env
    ),
    error = function(err) {
      cli_abort(
        c(
          paste(
            "An error occurred while calling the function {.fn {fun_name}} provided",
            "to the `get_terms_fun` argument."
          ),
          "This could be due to incorrect handling of input parameters inside {.fn {fun_name}}.",
          "Current arguments passed to {.fn {fun_name}}:",
          " - version: {version}",
          " - basket_select: {definition}",
          " - keep_id: {expect_grpid}",
          "Error message: {conditionMessage(err)}"
        )
      )
    }
  )

  assert_terms(
    terms,
    expect_grpname = expect_grpname,
    expect_grpid = expect_grpid,
    source_text = paste0(
      "object returned by calling get_terms_fun(basket_select = ",
      format(definition),
      ", version = ",
      dquote(version),
      ", keep_id = ",
      expect_grpid,
      ")"
    )
  )
  terms
}

#' Check required parameters for a basket
#'
#' If a basket (SMQ, SDG, ....) are requested, the version and a function to access the
#' database must be provided. The function checks these requirements.
#'
#' @param version Version provided by user
#'
#' @param version_arg_name Name of the argument providing the version
#'
#' @param fun Function provided by user
#'
#' @param fun_arg_name Name of the argument providing the function
#'
#' @param queries Queries provide by user
#'
#' @param i Index of query being checked
#'
#' @return An error is issued if `version` or `fun` is null.
#'
#' @keywords internal
assert_db_requirements <- function(version, version_arg_name, fun, fun_arg_name, queries, i) {
  if (is.null(fun)) {
    msg <-
      c(
        "{.arg {fun_arg_name}} is not specified. This is expected for baskets.",
        "A basket is requested by query {.val {i}}:",
        capture.output(queries[[i]])
      )
    # set names for correct indentation
    names(msg) <- if_else(str_starts(msg, " "), " ", "")
    names(msg)[[2]] <- "i"
    cli_abort(msg)
  }
  if (is.null(version)) {
    msg <-
      c(
        "{.arg {version_arg_name}} is not specified. This is expected for baskets.",
        "A basket is requested by query {.val {i}}:",
        capture.output(queries[[i]])
      )
    # set names for correct indentation
    names(msg) <- if_else(str_starts(msg, " "), " ", "")
    names(msg)[[2]] <- "i"
    cli_abort(msg)
  }
}

#' Create an `query` object
#'
#' A `query` object defines a query, e.g., a Standard MedDRA Query (SMQ), a
#' Standardized Drug Grouping (SDG), or a customized query (CQ). It is used
#' as input to `create_query_data()`.
#'
#' @param prefix The value is used to populate `PREFIX` in the output
#'   dataset of `create_query_data()`, e.g., `"SMQ03"`
#'
#' @param name The value is used to populate `GRPNAME` in the output dataset
#'   of `create_query_data()`. If the `auto` keyword is specified, the variable
#'   is set to the name of the query in the SMQ/SDG database.
#'
#' @permitted A character scalar or the `auto` keyword. The `auto`
#'   keyword is permitted only for queries which are defined by an
#'   `basket_select()` object.
#'
#' @param id The value is used to populate `GRPID` in the output dataset of
#'   `create_query_data()`. If the `auto` keyword is specified, the variable is
#'   set to the id of the query in the SMQ/SDG database.
#'
#' @permitted A integer scalar or the `auto` keyword. The `auto`
#'   keyword is permitted only for queries which are defined by an
#'   `basket_select()` object.
#'
#' @param add_scope_num Determines if  `SCOPEN` in the output dataset
#'   of `create_query_data()` is populated
#'
#'   If the parameter is set to `TRUE`, the definition must be an `basket_select()`
#'   object.
#'
#'   *Default*: `FALSE`
#'
#' @permitted `TRUE`, `FALSE`
#'
#' @param definition Definition of terms belonging to the query
#'
#'   There are three different ways to define the terms:
#'
#'   * An `basket_select()` object is specified to select a query from the SMQ
#'     database.
#'
#'   * A data frame with columns `SRCVAR` and `TERMCHAR` or `TERMNUM` can
#'     be specified to define the terms of a customized query. The `SRCVAR`
#'     should be set to the name of the variable which should be used to select
#'     the terms, e.g., `"AEDECOD"` or `"AELLTCD"`. `SRCVAR` does not need
#'     to be constant within a query. For example a query can be based on
#'     `AEDECOD` and `AELLT`.
#'
#'     If `SRCVAR` refers to a character variable, `TERMCHAR` should be set
#'     to the value the variable. If it refers to a numeric variable, `TERMNUM`
#'     should be set to the value of the variable. If only character variables
#'     or only numeric variables are used, `TERMNUM` or `TERMCHAR` respectively
#'     can be omitted.
#'
#'   * A list of data frames and `basket_select()` objects can be specified to
#'   define a customized query based on custom terms and SMQs. The data frames
#'   must have the same structure as described for the previous item.
#'
#' @permitted an `basket_select()` object, a
#'   data frame, or a list of data frames and `basket_select()` objects.
#'
#'
#' @seealso [create_query_data()], [basket_select()],
#'   `vignette("queries_dataset")`
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class `query`.
#'
#' @examples
#' # create a query for an SMQ
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' # create a query for a SMQ
#' query(
#'   prefix = "SMQ02",
#'   id = auto,
#'   definition = basket_select(
#'     name = "Pregnancy and neonatal topics (SMQ)",
#'     scope = "NARROW",
#'     type = "smq"
#'   )
#' )
#'
#' # create a query for an SDG
#' query(
#'   prefix = "SDG01",
#'   id = auto,
#'   definition = basket_select(
#'     name = "5-aminosalicylates for ulcerative colitis",
#'     scope = NA_character_,
#'     type = "sdg"
#'   )
#' )
#'
#' # creating a query for a customized query
#' cqterms <- tribble(
#'   ~TERMCHAR, ~TERMNUM,
#'   "APPLICATION SITE ERYTHEMA", 10003041L,
#'   "APPLICATION SITE PRURITUS", 10003053L
#' ) %>%
#'   mutate(SRCVAR = "AEDECOD")
#'
#' query(
#'   prefix = "CQ01",
#'   name = "Application Site Issues",
#'   definition = cqterms
#' )
#'
#' # creating a customized query based on SMQs and additional terms
#' query(
#'   prefix = "CQ03",
#'   name = "Special issues of interest",
#'   definition = list(
#'     cqterms,
#'     basket_select(
#'       name = "Pregnancy and neonatal topics (SMQ)",
#'       scope = "NARROW",
#'       type = "smq"
#'     ),
#'     basket_select(
#'       id = 8050L,
#'       scope = "BROAD",
#'       type = "smq"
#'     )
#'   )
#' )
query <- function(prefix,
                  name = auto,
                  id = NULL,
                  add_scope_num = FALSE,
                  definition = NULL) {
  out <- list(
    prefix = prefix,
    name = enexpr(name),
    id = enexpr(id),
    add_scope_num = add_scope_num,
    definition = definition
  )
  # evaluate to ensure that name contains the quoted symbol auto or a character
  # string
  if (!is_auto(out$name) && !is_missing(out$name)) {
    out$name <- eval_tidy(out$name)
  }
  # evaluate to ensure that id contains the quoted symbol auto or a number
  if (!is_auto(out$id)) {
    out$id <- eval_tidy(out$id)
  }
  class(out) <- c("query", "source", "list")
  validate_query(out)
}

#' Validate an object is indeed a `query` object
#'
#' @param obj An object to be validated.
#'
#' @seealso [query()]
#'
#' @return The original object.
#' @keywords internal
validate_query <- function(obj) {
  assert_s3_class(obj, "query")
  values <- unclass(obj)
  prefix <- values$prefix
  assert_character_scalar(prefix)

  if (!is_auto(values$name)) {
    name <- values$name
    assert_character_scalar(name)
  }

  if (!is_auto(values$id)) {
    id <- values$id
    assert_integer_scalar(id, optional = TRUE)
  }

  scope <- values$scope
  assert_character_scalar(scope,
    values = c("BROAD", "NARROW", NA_character_),
    optional = TRUE
  )

  add_scope_num <- values$add_scope_num
  assert_logical_scalar(add_scope_num,
    optional = TRUE
  )

  if (inherits(values$definition, "basket_select")) {
    validate_basket_select(values$definition)
  } else if (is.data.frame(values$definition) || is.list(values$definition)) {
    if (is_auto(values$name)) {
      cli_abort(
        c(
          "The auto keyword can be used for baskets only.",
          i = "It was provided for the {.var name} element."
        )
      )
    }
    if (is_auto(values$id)) {
      cli_abort(
        c(
          "The auto keyword can be used for baskets only.",
          i = "It was provided for the {.var id} element."
        )
      )
    }
    if (is.data.frame(values$definition)) {
      assert_terms(
        terms = values$definition,
        source_text = "the data frame provided for the `definition` element"
      )
    } else {
      is_valid <-
        map_lgl(values$definition, is.data.frame) |
          map_lgl(values$definition, inherits, "basket_select")
      if (!all(is_valid)) {
        info_msg <- paste(
          sprintf(
            "\u2716 Element %s is %s",
            which(!is_valid),
            map_chr(values$definition[!is_valid], what_is_it)
          ),
          collapse = "\n"
        )
        err_msg <- c(
          paste(
            "Each element of the list in the {.var definition} field must be a data frame",
            "or an object of class {.cls basket_select} but the following are not:"
          ),
          sprintf(
            "Element %s is {.obj_type_friendly {values$definition[[%s]]}}.",
            which(!is_valid),
            which(!is_valid)
          )
        )
        cli_abort(err_msg)
      }

      for (i in seq_along(values$definition)) {
        if (is.data.frame(values$definition[[i]])) {
          assert_terms(
            terms = values$definition[[i]],
            source_text = paste0("the ", i, "th element of the definition field")
          )
        }
      }
    }
  } else {
    cli_abort(
      c(
        paste(
          "{.var definition} expects a {.cls basket_select} object, a data frame,",
          "or a list of data frames and {.cls basket_select} objects."
        ),
        i = "An object of the following class was provided: {.cls {class(values$definition)}}"
      )
    )
  }
  obj
}

#' Asserts Requirements for Terms for Queries
#'
#' The function checks the requirements for terms for queries provided by the
#' user. The terms could have been provided directly in the query definition or
#' via a user provided function for accessing a SMQ or SDG database.
#'
#' @param terms Terms provided by user
#'
#' @param expect_grpname Is the `GRPNAME` column expected?
#'
#' @param expect_grpid Is the `GRPID` column expected?
#'
#' @param source_text Text describing the source of the terms, e.g., `"the data
#'   frame provided for the `definition` element"`.
#'
#' @return An error is issued if
#'
#' - `terms` is not a data frame,
#' - `terms` has zero observations,
#' - the `SRCVAR` variable is not in `terms`,
#' - neither the `TERMCHAR` nor the `TERMNUM` variable is in `terms`,
#' - `expect_grpname == TRUE` and the `GRPNAME` variable is not in `terms`,
#' - `expect_grpid == TRUE` and the `GRPID` variable is not in `terms`,
#'
#' @examples
#'
#' try(
#'   assert_terms(
#'     terms = 42,
#'     source_text = "object provided by the `definition` element"
#'   )
#' )
#'
#' @keywords internal
assert_terms <- function(terms,
                         expect_grpname = FALSE,
                         expect_grpid = FALSE,
                         source_text) {
  assert_data_frame(
    terms,
    message = paste0(
      source_text,
      " is not a data frame but {.obj_type_friendly {terms}}."
    )
  )

  if (nrow(terms) == 0) {
    cli_abort(paste0(
      source_text,
      " does not contain any observations."
    ))
  }

  required_vars <- exprs(SRCVAR)
  if (expect_grpname) {
    required_vars <- exprs(!!!required_vars, GRPNAME)
  }
  if (expect_grpid) {
    required_vars <- exprs(!!!required_vars, GRPID)
  }
  assert_data_frame(
    terms,
    required_vars = required_vars,
    message = paste0(
      "Required variable{?s} {.var {missing_vars}} {?is/are} missing in ",
      source_text,
      "."
    )
  )
  vars <- names(terms)
  if (!"TERMCHAR" %in% vars && !"TERMNUM" %in% vars) {
    cli_abort(
      c(
        "Variable {.var TERMCHAR} or {.var TERMNUM} is required.",
        paste0("None of them is in ", source_text, "."),
        i = "Provided variables: {.var {vars}}"
      )
    )
  }
}

#' Create a `basket_select` object
#'
#' @param name Name of the query used to select the definition of the query from
#'   the company database.
#'
#' @param id Identifier of the query used to select the definition of the query
#'   from the company database.
#'
#' @param scope Scope of the query used to select the definition of the query
#'   from the company database.
#'
#' @permitted `"BROAD"`, `"NARROW"`, `NA_character_`
#'
#' @param type The type argument expects a character scalar. It is passed to the
#' company specific get_terms() function such that the function can determine
#' which sort of basket is requested
#'
#' @param ... Any number of *named* function arguments. Can be used to pass in company
#' specific conditions or flags that will then be used in user-defined function that is
#' passed into argument `get_terms_fun` for function `create_query_data()`.
#'
#' @details Exactly one of `name` or `id` must be specified.
#'
#' @return An object of class `basket_select`.
#'
#'
#' @seealso [create_query_data()], [query()]
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @export
basket_select <- function(name = NULL,
                          id = NULL,
                          scope = NULL,
                          type,
                          ...) {
  args <- eval(substitute(alist(...)))
  if (length(args) == 0L) {
    out <- list(
      name = name,
      id = id,
      scope = scope,
      type = type
    )
  } else {
    if (!is_named(args)) {
      cli_abort("All arguments inside {.arg ...} must be named")
    }
    out <- list(
      name = name,
      id = id,
      scope = scope,
      type = type,
      ...
    )
  }

  class(out) <- c("basket_select", "source", "list")
  validate_basket_select(out)
}

#' Validate an object is indeed a `basket_select` object
#'
#' @param obj An object to be validated.
#'
#' @seealso [basket_select()]
#'
#' @keywords internal
#'
#' @return The original object.
validate_basket_select <- function(obj) {
  assert_s3_class(obj, "basket_select")
  values <- unclass(obj)
  name <- values$name
  assert_character_scalar(name,
    optional = TRUE
  )
  id <- values$id
  assert_integer_scalar(id,
    optional = TRUE
  )
  scope <- values$scope
  assert_character_scalar(scope,
    values = c("BROAD", "NARROW", NA_character_)
  )

  if (is.null(values$id) && is.null(values$name)) {
    cli_abort("Either {.var id} or {.var name} has to be non null.")
  }
  if (!is.null(values$id) && !is.null(values$name)) {
    cli_abort("Either {.var id} or {.var name} has to be null.")
  }
  obj
}

#' Returns a Character Representation of a `basket_select()` Object
#'
#' The function returns a character representation of a `basket_select()` object.
#' It can be used for error messages for example.
#'
#' @param x A `basket_select()` object
#'
#' @param ... Not used
#'
#' @return A character representation of the `basket_select()` object
#'
#'
#' @seealso [basket_select()]
#'
#' @keywords internal
#' @family internal
#'
#' @export
#'
#' @examples
#'
#' format(basket_select(id = 42, scope = "NARROW", type = "smq"))
format.basket_select <- function(x, ...) {
  all_arg_names <- names(x)

  formvar <- list()

  for (i in seq_len(length(all_arg_names))) {
    is_numeric_class <- map_lgl(x[i], inherits, "numeric") | map_chr(x[i], typeof) == "numeric"

    if (is.character(x[[i]]) && length(x[[i]]) <= 1 && !is.na(x[[i]])) {
      formvar[i] <- paste(all_arg_names[i], "=", dquote(x[[i]]))
    } else if (length(x[[i]]) <= 1 || typeof(x[[i]]) == "language") {
      formvar[i] <- paste(
        all_arg_names[i],
        "=",
        format(x[[i]]) %>% trimws() %>% paste(collapse = " ")
      )
    } else {
      if (typeof(x[[i]]) == "list") {
        obj_name <- "list"
      } else {
        obj_name <- "c"
      }
      formvar[i] <- paste0(all_arg_names[i], " = ", obj_name, "(...)")
    }
  }

  allvars <- paste(formvar, collapse = ", ")

  paste0(
    "basket_select(",
    allvars,
    ")"
  )
}
