#' Creates queries dataset from SMQ and SDG database or custom definitions as
#' input dataset for `derive_vars_query()`
#'
#' Creates queries dataset from definitions in SMQ and SDG database or custom
#' definitions as input dataset for `derive_vars_query()`.
#'
#' @param queries List of queries
#'
#'   A list of `query()` objects is expected.
#'
#' @param meddra_version MedDRA version
#'
#'   The MedDRA version used for coding the terms in the AE dataset should be
#'   specified. If any of the queries is a SMQ or a customized query including a
#'   SMQ, the parameter needs to be specified.
#'
#'   *Permitted Values*: A character string
#'
#' @param whodd_version WHO Drug Dictionary version
#'
#'   The version of the WHO Drug Dictionary used for coding the terms in the CM
#'   dataset should be specified. If any of the queries is a SDG, the parameter
#'   needs to be specified.
#'
#'   *Permitted Values*: A character string
#'
#' @param get_smq_fun Function which returns the terms of a SMQ
#'
#'   For each query specified for the `queries` parameter which refers to a SMQ
#'   the specified function is called to retrieve the terms defining the query.
#'
#'   The function must return a dataset with all the terms defining the SMQ. The
#'   output dataset must contain the following variables.
#'
#'   - `TERM_LEVEL`
#'   - `TERM_NAME` if the variable `TERM_LEVEL` is referring to is character
#'   - `TERM_ID` if the variable `TERM_LEVEL` is referring to is numeric
#'   - `QUERY_NAME`
#'
#'   The function must provide the following parameters
#'
#'   - `smq_select`: A `smq_select()` object.
#'   - `version`: The MedDRA version. The value specified for the
#'   `meddra_version` in the `create_query_data()` call is passed to this
#'   parameter.
#'   - `keep_id`: If set to `TRUE`, the output dataset must contain the
#'   `QUERY_ID` variable
#'   - `temp_env`: A temporary environment is passed to this parameter. It can
#'   be used to store data which is used for all SMQs in the
#'   `create_query_data()` call. For example if the SMQs need to be read from a
#'   database all SMQs can be read and stored in the environment when the first
#'   SMQ is handled. For the other SMQs the terms can be retrieved from the
#'   environment instead of accessing the database again.
#'
#' @param get_sdg_fun Function which returns the terms of a SDG
#'
#'   For each query specified for the `queries` parameter which refers to a SDG
#'   the specified function is called to retrieve the terms defining the query.
#'
#'   The function must return a dataset with all the terms defining the SDG. The
#'   output dataset must contain the following variables.
#'
#'   - `TERM_LEVEL`
#'   - `TERM_NAME` if the variable `TERM_LEVEL` is referring to is character
#'   - `TERM_ID` if the variable `TERM_LEVEL` is referring to is numeric
#'   - `QUERY_NAME`
#'
#'   The function must provide the following parameters
#'
#'   - `sdg_select`: A `sdg_select()` object.
#'   - `version`: The WHO drug dictionary version. The value specified for the
#'   `whodd_version` in the `create_query_data()` call is passed to this
#'   parameter.
#'   - `keep_id`: If set to `TRUE`, the output dataset must contain the
#'   `QUERY_ID` variable
#'   - `temp_env`: A temporary environment is passed to this parameter. It can
#'   be used to store data which is used for all SDGs in the
#'   `create_query_data()` call. For example if the SDGs need to be read from a
#'   database all SDGs can be read and stored in the environment when the first
#'   SDG is handled. For the other SDGs the terms can be retrieved from the
#'   environment instead of accessing the database again.
#'
#' @details For each query the terms belonging to the query are determined with
#'   respect to the `definition` element.
#'   * If it is a `smq_select()` object, the terms are read from the SMQ
#'   database by calling the function specified for the `get_smq_fun` parameter.
#'   * If it is a `sdg_select()` object, the terms are read from the SDG
#'   database by calling the function specified for the `get_sdg_fun` parameter.
#'   * If it is a data frame, the terms stored in the data frame are used.
#'   * If it is a list of data frames and `smq_select()` objects, all terms from
#'   the dataframes and all terms read from the SMQ database referrenced by the
#'   `smq_select` objects are collated.
#'
#'   The following variables are created:
#'
#'   * `VAR_PREFIX`: Prefix of the variables to be created by
#'   `derive_vars_query()` as specified by the `prefix` element.
#'   * `QUERY_NAME`: Name of the query as specified by the `name` element.
#'   * `QUERY_ID`: Id of the query as specified by the `id` element. If the `id`
#'   element is not specified for a query, the variable is set to `NA`. If the
#'   `id` element is not specified for any query, the variable is not created.
#'   * `QUERY_SCOPE`: scope of the query as specified by the `scope` element of
#'   the `smq_select()` object. For queries not defined by a `smq_select()`
#'   object, the variable is set to `NA`. If none of the queries is defined by a
#'   `smq_select()` object, the variable is not created.
#'   * `QUERY_SCOPE_NUM`: numeric scope of the query. It is set to `1` if the
#'   scope is broad. Otherwise it is set to '2'. If the `add_scope_num` element
#'   equals `FALSE`, the variable is set to `NA`. If the `add_scope_num` element
#'   equals `FALSE` of none of the queries is a SMQ , the variable is not
#'   created.
#'   * `TERM_LEVEL`: Name of the variable used to identify the terms.
#'   * `TERM_NAME`: Value of the term variable if it is a character variable.
#'   * `TERM_ID`: Value of the term variable if it is a numeric variable.
#'
#' @author Stefan Bundfuss
#'
#' @return A dataset to be used as input for `derive_vars_query()`
#'
#' @keywords adae adcm user_utility
#'
#' @seealso [derive_vars_query()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(magrittr)
#' library(dplyr)
#'
#' # creating a query dataset for a customized query
#' cqterms <- tribble(
#'   ~TERM_NAME, ~TERM_ID,
#'   "APPLICATION SITE ERYTHEMA", 10003041L,
#'   "APPLICATION SITE PRURITUS", 10003053L) %>%
#'   mutate(TERM_LEVEL = "AEDECOD")
#'
#' cq <- query(prefix = "CQ01",
#'             name = "Application Site Issues",
#'             definition = cqterms)
#'
#' create_query_data(queries = list(cq))
#'
#' # create a query dataset for SMQs
#' pregsmq <- query(
#'   prefix = "SMQ02",
#'   id = auto,
#'   definition = smq_select(name = "Pregnancy and neonatal topics (SMQ)",
#'                           scope = "NARROW"))
#'
#' pneusmq <- query(prefix = "SMQ04",
#'                  definition = smq_select(id = 8050L,
#'                                          scope = "BROAD"))
#'
#' \dontrun{
#' create_query_data(queries = list(pregsmq, pneusmq),
#'                   meddra_version = "20.1")
#' }
#'
#' # create a query dataset for SDGs
#' sdg <- query(
#'   prefix = "SDG01",
#'   id = auto,
#'   definition = sdg_select(
#'     name = "5-aminosalicylates for ulcerative colitis"
#'   )
#' )
#' \dontrun{
#' create_query_data(queries = list(sdg),
#'                   whodd_version = "2019_09")
#' }
#'
#' # creating a query dataset for a customized query including SMQs
#' \dontrun{
#' create_query_data(queries = list(pregsmq, cqterms),
#'                   meddra_version = "20.1")
#' }

create_query_data <- function(queries,
                              meddra_version = NULL,
                              whodd_version = NULL,
                              get_smq_fun = NULL,
                              get_sdg_fun = NULL) {
  # check parameters
  assert_character_scalar(meddra_version, optional = TRUE)
  assert_character_scalar(whodd_version, optional = TRUE)
  assert_function(get_smq_fun,
                  params = c("smq_select", "version", "keep_id", "temp_env"),
                  optional = TRUE)
  assert_function(get_sdg_fun,
                  params = c("sdg_select", "version", "keep_id", "temp_env"),
                  optional = TRUE)

  walk(queries, validate_query)

  # read queries
  temp_env <- new.env(parent = emptyenv())
  query_data <- vector("list", length(queries))
  for (i in seq_along(queries)) {
    # get term names and term variable
    if (inherits(queries[[i]]$definition, "smq_select")) {
      # query is a SMQ
      assert_db_requirements(version = meddra_version,
                             fun = get_smq_fun,
                             queries = queries,
                             i = i,
                             type = "SMQ")

      query_data[[i]] <- call_user_fun(
        get_smq_fun(
          smq_select = queries[[i]]$definition,
          version = meddra_version,
          keep_id = !is.null(queries[[i]]$id),
          temp_env = temp_env
        )
      )
      assert_terms(
        query_data[[i]],
        expect_query_name = TRUE,
        expect_query_id = !is.null(queries[[i]]$id),
        source_text = paste0(
          "object returned by calling get_smq_fun(smq_select = ",
          format(queries[[i]]$definition),
          ", version = ",
          dquote(meddra_version),
          ", keep_id = ",!is.null(queries[[i]]$id),
          ")"
        )
      )
      query_data[[i]] <- mutate(query_data[[i]],
                                QUERY_SCOPE = queries[[i]]$definition$scope)
      if (queries[[i]]$add_scope_num) {
        query_data[[i]] <-
          mutate(query_data[[i]],
                 QUERY_SCOPE_NUM = if_else(QUERY_SCOPE == "BROAD", 1, 2))
      }
    }
    else if (inherits(queries[[i]]$definition, "sdg_select")) {
      # query is a SDG
      assert_db_requirements(version = whodd_version,
                             fun = get_sdg_fun,
                             queries = queries,
                             i = i,
                             type = "SDG")
      query_data[[i]] <- call_user_fun(
        get_sdg_fun(
          sdg_select = queries[[i]]$definition,
          version = whodd_version,
          keep_id = !is.null(queries[[i]]$id),
          temp_env = temp_env
        )
      )
      assert_terms(
        query_data[[i]],
        expect_query_name = TRUE,
        expect_query_id = !is.null(queries[[i]]$id),
        source_text = paste0(
          "object returned by calling get_sdg_fun(sdg_select = ",
          format(queries[[i]]$definition),
          ", version = ",
          dquote(whodd_version),
          ", keep_id = ",!is.null(queries[[i]]$id),
          ")"
        )
      )
    }
    else if (is.data.frame(queries[[i]]$definition)) {
      # query is a customized query
      query_data[[i]] <- queries[[i]]$definition
    }
    else if (is.list(queries[[i]]$definition)) {
      # query is defined by customized queries and SMQs
      definition <- queries[[i]]$definition
      terms <- vector("list", length(definition))
      for (j in seq_along(definition)) {
        if (is.data.frame(definition[[j]])) {
          terms[[j]] <- definition[[j]]
        }
        else {
          assert_db_requirements(
            version = meddra_version,
            fun = get_smq_fun,
            queries = queries,
            i = i,
            type = "SMQ"
          )
          terms[[j]] <-
            call_user_fun(get_smq_fun(
              smq_select = definition[[j]],
              version = meddra_version,
              temp_env = temp_env
            ))
          assert_terms(
            terms[[j]],
            source_text = paste0(
              "object returned by calling get_smq_fun(smq_select = ",
              format(definition[[j]]),
              ", version = ",
              dquote(meddra_version),
              ")"
            )
          )

        }
      }
      query_data[[i]] <- bind_rows(terms)
    }

    # add mandatory variables
    query_data[[i]] <- mutate(query_data[[i]],
                              VAR_PREFIX = queries[[i]]$prefix)

    if (!is_auto(queries[[i]]$name)) {
      query_data[[i]] <- mutate(query_data[[i]],
                                QUERY_NAME = queries[[i]]$name)
    }

    # add optional variables
    if (!is.null(queries[[i]]$id) && !is_auto(queries[[i]]$id)) {
      query_data[[i]] <- mutate(query_data[[i]],
                                QUERY_ID = queries[[i]]$id)
    }
  }
  bind_rows(query_data)
}

assert_db_requirements <- function(version, fun, queries, i, type) {
  if (is.null(fun)) {
    msg <-
      paste0(
        arg_name(substitute(fun)),
        " is not specified. This is expected for ",
        type,
        "s.\n",
        "A ", type, " is requested by query ",
        i,
        ":\n",
        paste(capture.output(str(queries[[i]])),
              collapse = "\n")
      )
    abort(msg)
  }
  if (is.null(version)) {
    msg <-
      paste0(
        arg_name(substitute(version)),
        " is not specified. This is expected for ",
        type,
        "s.\n",
        "A ",
        type,
        " is requested by query ",
        i,
        ":\n",
        paste(capture.output(str(queries[[i]])),
              collapse = "\n")
      )
    abort(msg)
  }
}

#' Create an `query` object
#'
#' An `query` object defines a query, e.g., a Standard MedDRA Query (SMQ), a
#' Standardised Drug Grouping (SDG), or a customized query (CQ).
#'
#' @param prefix Prefix of the new variables, e.g., `"SMQ03"`
#'
#' @param name The variable `<prefix>NAM` is set to the specified value for
#'   observations matching the query definition. If the `auto` keyword is
#'   specified, the variable is set to the name of the query in the SMQ/SDG
#'   database.
#'
#'   *Permitted Values*: A character scalar or the `auto` keyword. The `auto`
#'   keyword is permitted only for queries which are defined by a `smq_select()`
#'   or `sdg_select()` object.
#'
#' @param id The variable `<prefix>CD` is set to the specified value for
#'   observations matching the query definition if the parameter is specified.
#'   Otherwise the variable is not created. If the `auto` keyword is specified,
#'   the variable is set to the id of the query in the SMQ/SDG database.
#'
#'   *Permitted Values*: A integer scalar or the `auto` keyword. The `auto`
#'   keyword is permitted only for queries which are defined by a `smq_select()`
#'   or `sdg_select()` object.
#'
#'  @param add_scope_num The variable `<prefix>SCN` is added to the output
#'   dataset if the parameter is specified. Otherwise the variable is not
#'   created.
#'
#'   The value is set to `1` if the scope is broad. Otherwise, it is set to `2`.
#'
#'   If the parameter is set to `TRUE`, the definition must be a `smq_select()`
#'   object.
#'
#'   *Default*: `FALSE`
#'   *Permitted Values*: `TRUE`, `FALSE`
#'
#' @param definition Definition of terms belonging to the query
#'
#'   There are three different ways to define the terms:
#'
#'   * A `smq_select()` object is specified to select a query from the SMQ
#'     database. The database contains Standard MedDRA Queries (SMQs) and AE
#'     Grouping Terms (AEGTs).
#'
#'   * A `sdg_select()` object is specified to select a query from the SDG
#'     database.
#'
#'   * A data frame with columns `TERM_LEVEL` and `TERM_NAME` or `TERM_ID` can
#'     be specified to define the terms of the query. The `TERM_LEVEL` should be
#'     set to the name of the variable which should be used to select the terms,
#'     e.g., `"AEDECOD"` or `"AELLTCD"`. `TERM_LEVEL` does not need to be
#'     constant within a query. For example a query can be based on `AEDECOD`
#'     and `AELLT`.
#'
#'     If `TERM_LEVEL` refers to a character variable, `TERM_NAME` should be set
#'     to the value the variable. If it refers to a numeric variable, `TERM_ID`
#'     should be set to the value of the variable. If only character variables
#'     or only numeric variables are used, `TERM_ID` or `TERM_NAME` respectively
#'     can be omitted.
#'
#' @details
#'
#' The definition includes
#'
#' * the variables which should be created for the query, e.g., `SMQ02NAM` and
#' `SMQ02SC`,
#'
#' * the values which should be assigned to the variables, and
#'
#' * which terms belong to the query.
#'
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "query".
query <- function(prefix,
                  name = auto,
                  id = NULL,
                  add_scope_num = FALSE,
                  definition = NULL) {
  out <- list(
    prefix = prefix,
    name = enquo(name),
    id = enquo(id),
    add_scope_num = add_scope_num,
    definition = definition
  )
  # evaluate to ensure that name contains the quoted symbol auto or a character
  # string
  if (!is_auto(out$name) && !quo_is_missing(out$name)) {
    out$name <- eval_tidy(out$name)
  }
  # evaluate to ensure that id contains the quoted symbol auto or a number
  if (!is_auto(out$id)) {
    out$id <- eval_tidy(out$id)
  }
  class(out) <- c("query", "list")
  validate_query(out)
}

#' Validate an object is indeed a `query` object
#'
#' @param obj An object to be validated.
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @return The original object.
validate_query <- function(obj) {
  assert_that(inherits(obj, "query"))
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
                          values = c("BROAD", "NARROW"),
                          optional = TRUE)

  add_scope_num <- values$add_scope_num
  assert_logical_scalar(add_scope_num,
                        optional = TRUE)
  if (add_scope_num && !inherits(values$definition, "smq_select")) {
    abort("`add_scope_num == TRUE` must be used for SMQs only.")
  }

  if (inherits(values$definition, "smq_select")) {
    validate_smq_select(values$definition)
  }
  else if (inherits(values$definition, "sdg_select")) {
    validate_sdg_select(values$definition)
  }
  else if (is.data.frame(values$definition) || is.list(values$definition)) {
    if (is_auto(values$name)) {
      abort(
        paste0(
          "The auto keyword can be used for SMQs and SDGs only.\n",
          "It was provided for the name element."
        )
      )
    }
    if (is_auto(values$id)) {
      abort(
        paste0("The auto keyword can be used for SMQs and SDGs only.\n",
               "It was provided for the id element."))
    }
    if (is.data.frame(values$definition)) {
      assert_terms(terms = values$definition,
                   source_text = "the data frame provided for the `definition` element")
    }
    else {
      is_valid <- map_lgl(values$definition, is.data.frame) | map_lgl(values$definition, inherits, "smq_select")
      if (!all(is_valid)) {
        info_msg <- paste(
          sprintf("\u2716 Element %s is %s", which(!is_valid), map_chr(values$definition[!is_valid], what_is_it)),
          collapse = "\n"
        )
        err_msg <- sprintf(
          "Each element of the list in the definition field must be a data frame or an object of class `smq_select` but the following are not:\n%s",
          info_msg
        )
        abort(err_msg)
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
  }
  else {
    abort(
      paste0(
        "`definition` expects a `smq_select` or `sdg_select` object, a data frame, or a list of data frames and `smq_select` objects.\n",
        "An object of the following class was provided: ",
        class(values$definition)
      )
    )
  }
  obj
}

assert_terms <- function(terms,
                         expect_query_name = FALSE,
                         expect_query_id = FALSE,
                         source_text) {
  if (!is.data.frame(terms)) {
    abort(paste0(source_text,
                 " is not a data frame but ",
                 what_is_it(terms),
                 "."))
  }

  if (nrow(terms) == 0) {
    abort(paste0(
          source_text,
          " does not contain any observations."))
  }

  vars <- names(terms)
  if (!"TERM_LEVEL" %in% vars) {
    abort(
      paste0("Required variable `TERM_LEVEL` is missing in ",
            source_text,
            ".")
    )
  }
  if (expect_query_name) {
    if (!"QUERY_NAME" %in% vars) {
      abort(
        paste0("Required variable `QUERY_NAME` is missing in ",
               source_text,
               ".")
      )
    }
  }
  if (expect_query_id) {
    if (!"QUERY_ID" %in% vars) {
      abort(
        paste0("Required variable `QUERY_ID` is missing in ",
               source_text,
               ".")
      )
    }
  }
  if (!"TERM_NAME" %in% vars & !"TERM_ID" %in% vars) {
    abort(
      paste0(
        "Variable `TERM_NAME` or `TERM_ID` is required.\n",
        "None of them is in ",
        source_text,
        ".\n",
        "Provided variables: ",
        enumerate(vars)
      )
    )
  }
}

#' Create an `smq_select` object
#'
#' @param name Name of the query used to select the definition of the query.
#'
#' @param id Identifier of the query used to select the definition of the query.
#'
#' @param scope Scope of the query used to select the definition of the query.
#'
#'   *Permitted Values*: `"BROAD"`, `"NARROW"`
#'
#' @details Exactly one of `name` or `id` must be specified.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "smq_select".
smq_select <- function(name = NULL,
                       id = NULL,
                       scope = NULL) {
  out <- list(
    name = name,
    id = id,
    scope = scope
  )
  class(out) <- c("smq_select", "list")
  validate_smq_select(out)
}

#' Validate an object is indeed a `smq_select` object
#'
#' @param obj An object to be validated.
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @return The original object.
validate_smq_select <- function(obj) {
  assert_that(inherits(obj, "smq_select"))
  values <- unclass(obj)
  name <- values$name
  assert_character_scalar(name,
                          optional = TRUE)
  id <- values$id
  assert_integer_scalar(id,
                        optional = TRUE)
  scope <- values$scope
  assert_character_scalar(scope,
                          values = c("BROAD", "NARROW"))

  if (is.null(values$id) && is.null(values$name)) {
    abort("Either id or name has to be non null.")
  }
  if (!is.null(values$id) && !is.null(values$name)) {
    abort("Either id or name has to be null.")
  }
  obj
}

format.smq_select <- function(x, ...) {
  paste0("smq_select(name = ",
         dquote(x$name),
         ", id = ",
         format(x$id),
         ", scope = ",
         dquote(x$scope),
         ")")
}

#' Create an `sdg_select` object
#'
#' @param name Name of the query used to select the definition of the query
#'   from the central database.
#'
#'   The query with `GROUP_NAME == "<specified value>"` is selected.
#'
#' @param id Identifier of the query used to select the definition of the query
#'   from the central database.
#'
#'   The query with `GROUP_ID == <specified value>` is selected.
#'
#' @details Exactly one `name` or `id` must be specified.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "sdg_select".
sdg_select <- function(name=NULL,
                       id = NULL) {
  out <- list(
    name = name,
    id = id
  )
  class(out) <- c("sdg_select", "list")
  validate_sdg_select(out)
}

#' Validate an object is indeed a `sdg_select` object
#'
#' @param obj An object to be validated.
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @return The original object.
validate_sdg_select <- function(obj) {
  assert_that(inherits(obj, "sdg_select"))
  values <- unclass(obj)
  name <- values$name
  assert_character_scalar(name,
                          optional = TRUE)
  id <- values$id
  assert_integer_scalar(id,
                        optional = TRUE)
  if (is.null(values$id) && is.null(values$name)) {
    abort("Either id or name has to be non null.")
  }
  if (!is.null(values$id) && !is.null(values$name)) {
    abort("Either id or name has to be null.")
  }
  obj
}

format.sdg_select <- function(x, ...) {
  paste0("sdg_select(name = ",
         dquote(x$name),
         ", id = ",
         format(x$id),
         ")")
}
