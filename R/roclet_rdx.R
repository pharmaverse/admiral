#' Roclet Extending the Standard rd Roclet
#'
#' This roclet extends the standard `rd` roclet by allowing
#' - to add permitted values and default values to the `@param` tag and
#' - to add a caption and a description to examples.
#'
#' The following tags are supported:
#'
#' - `@permitted`: Permitted values for the argument. Permitted value description
#'     which are used for several arguments/functions can be stored in
#'     `inst/roxygen/rdx_meta.R`. For example:
#'     ```
#'     list(
#'       rdx_permitted_values = list(
#'         mode = "`\"first\"`, `\"last\"`",
#'         msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`"
#'       )
#'     )
#'     ```
#'     The reference to the permitted values is done by specifying the name of
#'     the list element in square brackets, e.g., `@permitted [mode]`.
#' - `@default`: Default value for the argument. By default the default value
#'   from the function formals is displayed. This can be overwritten by using
#'   the `@default` tag.
#'
#' - `@examplesx`: This tag can be used to mark the beginning of the examples
#'    section but doesn't affect the output, i.e., it can be omitted.
#'
#' - `@caption`: Caption for the example. The caption is displayed as a subsection
#'   in the examples section. The caption can be followed by an arbitrary number
#'   of `@info` and `@code` tags.
#' - `@info`: Description of the example.
#' - `@code`: Code of the example.
#'
#' To use the roclet call `roxygen2::roxygenise(roclets =
#' "admiral::rdx_roclet")` or add to the `DESCRIPTION` file:
#' ```
#' Roxygen: list(markdown = TRUE, roclets = c("collate", "namespace", "admiral::rdx_roclet"))
#' ```
#'
#' @keywords internal
#'
#' @export
rdx_roclet <- function() {
  out <- roclet("rdx")
  class(out) <- c("roclet_rdx", "roclet_rd", "roclet")
  out
}

#' @export
roclet_process.roclet_rdx <- function(x, blocks, env, base_path) {
  # read metadata for permitted values
  if (file.exists("./man/roxygen/rdx_meta.R")) {
    # nolint next
    rdx_permitted_values <- source("./man/roxygen/rdx_meta.R")$value$rdx_permitted_values
  } else {
    rdx_permitted_values <- NULL
  }
  # pre-process blocks to
  # - move @permitted and @default tags into @param tags and add defaults from function call
  # - transform @caption, @info and @code tags into @examplex tags
  blocks <- blocks %>%
    map(transform_param, rdx_permitted_values = rdx_permitted_values) %>%
    map(transform_examplesx)
  # call standard processing of rd roclet
  NextMethod()
}

transform_param <- function(block, rdx_permitted_values) {
  tags <- block$tags
  out_tags <- list()
  act_param <- list()
  if (length(block$call) >= 3) {
    defaults <- block$call[[3]][[2]]
  } else {
    defaults <- list()
  }
  for (i in seq_along(tags)) {
    if (tags[[i]]$tag == "param") {
      if (length(act_param) > 0) {
        out_tags <- c(out_tags, list(get_param_tag(act_param, defaults)))
      }
      act_param <- list(tag = tags[[i]], permitted = NULL, default = NULL)
    } else if (tags[[i]]$tag == "permitted") {
      if (!is.na(tags[[i]]$val$ref)) {
        ref_resolved <- rdx_permitted_values[[tags[[i]]$val$ref]]
        if (!is.null(ref_resolved)) {
          temp_tag <- tags[[i]]
          temp_tag$raw <- ref_resolved
          ref_resolved <- tag_markdown(temp_tag)$val
        }
        tags[[i]]$val$ref <- ref_resolved
      }
      act_param$permitted <- tags[[i]]$val %>%
        discard(is.na) %>%
        paste(collapse = " ")
    } else if (tags[[i]]$tag == "default") {
      act_param$default <- tags[[i]]$val
    } else {
      if (length(act_param) > 0) {
        out_tags <- c(out_tags, list(get_param_tag(act_param, defaults)))
        act_param <- list()
      }
      if (tags[[i]]$tag != "param") {
        out_tags <- c(out_tags, list(tags[[i]]))
      }
    }
  }
  block$tags <- out_tags
  block
}

get_param_tag <- function(act_param, defaults) {
  tag <- act_param$tag
  if (is.null(act_param$default)) {
    default_value <- defaults[[tag$val$name]]
    if (is_missing(default_value)) {
      default_value <- "none"
    } else {
      default_value <- paste0("\\code{", expr_deparse(default_value), "}")
    }
    act_param$default <- default_value
  }
  tag$val$description <- paste0(
    tag$val$description,
    "\n\n\\describe{\n",
    if_else(
      !is.null(act_param$permitted),
      paste0("\\item{Permitted values}{", act_param$permitted, "}\n"),
      ""
    ),
    "\\item{Default Value}{", act_param$default, "}\n}"
  )
  tag
}

transform_examplesx <- function(block) {
  tags <- block$tags
  out_tags <- list()
  act_example <- list()
  example_env <- new_environment(parent = global_env())
  for (i in seq_along(tags)) {
    if (tags[[i]]$tag == "caption") {
      if (length(act_example) > 0) {
        out_tags <- c(
          out_tags,
          list(roxy_tag("examplex", "generated", val = act_example))
        )
      }
      act_example <- list(caption = tags[[i]]$val, contents = NULL)
    } else if (tags[[i]]$tag == "info") {
      act_example$contents <- paste(act_example$contents, tags[[i]]$val, sep = "\n\n")
    } else if (tags[[i]]$tag == "code") {
      act_example$contents <- paste(
        act_example$contents,
        paste(
          "\\if{html}{\\out{<div class=\"sourceCode r\">}}\\preformatted{",
          execute_example(tags[[i]]$raw, env = example_env),
          "}\\if{html}{\\out{</div>}}",
          sep = ""
        ),
        sep = "\n\n"
      )
    } else {
      if (length(act_example) > 0) {
        out_tags <- c(
          out_tags,
          list(roxy_tag("examplex", "generated", val = act_example))
        )
        act_example <- list()
      }
      if (tags[[i]]$tag != "examplesx") {
        out_tags <- c(out_tags, list(tags[[i]]))
      }
    }
  }
  block$tags <- out_tags
  block
}

execute_example <- function(code, env = caller_env()) {
  expr_list <- parse(text = code)
  result <- NULL
  for (i in seq_along(expr_list)) {
    result <- c(result, as.character(attr(expr_list, "srcref")[[i]]))
    return_value <- withVisible(eval(expr_list[[i]], envir = env))
    if (return_value$visible) {
      result <- c(result, paste("#>", capture.output(print(return_value$value))))
    }
  }
  paste(result, collapse = "\n")
}

#' @export
roxy_tag_parse.roxy_tag_examplesx <- function(x) {
  x
}

#' @export
roxy_tag_parse.roxy_tag_examplex <- function(x) {
  tag_markdown(x)
}

#' @export
roxy_tag_rd.roxy_tag_examplex <- function(x, base_path, env) {
  rd_section("examplex", x$val)
}

#' @export
merge.rd_section_examplex <- function(x, y, ...) {
  stopifnot(identical(class(x), class(y)))
  rd_section(
    x$type,
    list(
      caption = c(x$value$caption, y$value$caption),
      contents = c(x$value$contents, y$value$contents)
    )
  )
}

#' @export
format.rd_section_examplex <- function(x, ...) {
  paste0(
    "\\section{Examples}{\n",
    paste0("\\subsection{", x$value$caption, "}{", x$value$contents, "}", collapse = "\n"),
    "}\n"
  )
}

#' @export
roxy_tag_parse.roxy_tag_caption <- function(x) {
  tag_markdown(x)
}

#' @export
roxy_tag_parse.roxy_tag_info <- function(x) {
  tag_markdown(x)
}

#' @export
roxy_tag_parse.roxy_tag_code <- function(x) {
  x
}
#' @export
roxy_tag_parse.roxy_tag_permitted <- function(x) {
  raw_parsed <- str_match(x$raw, "(?:\\[(.*)\\] *)?(.+)?")[, 2:3]
  x_text <- x
  if (is.na(raw_parsed[[2]])) {
    x_text$val <- raw_parsed[[2]]
  } else {
    x_text$raw <- raw_parsed[[2]]
    x_text <- tag_markdown(x_text)
  }
  x$val <- list(
    ref = raw_parsed[[1]],
    text = x_text$val
  )
  x
}

#' @export
roxy_tag_parse.roxy_tag_default <- function(x) {
  tag_markdown(x)
}
