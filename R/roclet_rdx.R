#' @export
rdx_roclet <- function() {
  out <- roclet("rdx")
  class(out) <- c("roclet_rdx", "roclet_rd", "roclet")
  out
}

#' @export
roclet_process.roclet_rdx <- function(x, blocks, env, base_path) {
  # print(blocks)
  # blocks[[1]]$tags[[2]]$tag <- "paramx"
  # print(blocks)
  blocks <- map(blocks, transform_examplesx)
  NextMethod()
  # Convert each block into a topic, indexed by filename
  # topics <- roxygen2:::RoxyTopics$new()
  #
  # for (block in blocks) {
  #   rd <- roxygen2:::block_to_rd(block, base_path, env)
  #   topics$add(rd, block)
  # }
  # roxygen2:::topics_process_family(topics, env)
  # roxygen2:::topics_process_inherit(topics, env)
  # # roxygen2:::topics$drop_invalid()
  # roxygen2:::topics_fix_params_order(topics)
  # roxygen2:::topics_add_default_description(topics)
  # roxygen2:::topics_add_package_alias(topics)
  #
  # topics$topics
}

transform_examplesx <- function(block) {
  tags <- block$tags
  out_tags <- list()
  act_example <- list()
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
          str_remove_all(tags[[i]]$raw, "(^\n+|\n+$)"),
          "}\\if{html}{\\out{</div>}}",
          sep = ""
        ),
        paste(
          "\\if{html}{\\out{<div class=\"sourceCode r\">}}\\Sexpr[stage=render,results=verbatim]{",
          str_remove_all(tags[[i]]$raw, "(^\n+|\n+$)"),
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

#' @export
roxy_tag_parse.roxy_tag_tip <- function(x) {
  tag_markdown(x)
}

#' @export
roxy_tag_rd.roxy_tag_tip <- function(x, base_path, env) {
  rd_section("tip", x$val)
}

#' @export
format.rd_section_tip <- function(x, ...) {
  paste0(
    "\\section{Tips and tricks}{\n",
    "\\itemize{\n",
    paste0("  \\item ", x$value, "\n", collapse = ""),
    "}\n",
    "}\n"
  )
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
  rd_section(x$type, list(caption = c(x$value$caption, y$value$caption), contents = c(x$value$contents, y$value$contents)))
}

#' @export
format.rd_section_examplex <- function(x, ...) {
  paste0(
    "\\section{Examples}{\n",
    paste0("\\subsection{", x$value$caption, "}{", x$value$contents, "}", collapse = "\n"),
    # paste0("\\subsection{X example}{", x$value, "}", collapse = ""),
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
