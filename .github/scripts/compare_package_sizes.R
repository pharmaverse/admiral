parse_args <- function(args) {
  defaults <- list(
    output_dir = file.path(getwd(), "reports"),
    package_name = "admiral",
    cran_mirror = "https://cran.r-project.org",
    github_repo = "pharmaverse/admiral",
    github_ref = "main"
  )

  if (length(args) == 0) {
    return(defaults)
  }

  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
  }

  parsed <- defaults

  for (index in seq(1, length(args), by = 2)) {
    key <- args[[index]]
    value <- args[[index + 1]]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument format: ", key, call. = FALSE)
    }

    parsed[[sub("^--", "", key)]] <- value
  }

  parsed
}

run_du <- function(args) {
  output <- system2(
    command = "du",
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste(c(sprintf("Command failed: du %s", paste(args, collapse = " ")), output), collapse = "\n"),
      call. = FALSE
    )
  }

  output[[1]]
}

measure_installed_package <- function(package_name, library_dir) {
  package_path <- find.package(package_name, lib.loc = library_dir)
  size_human <- sub("[[:space:]].*$", "", run_du(c("-sh", package_path)))
  size_bytes <- suppressWarnings(as.numeric(sub("[[:space:]].*$", "", run_du(c("-sb", package_path)))))

  if (is.na(size_bytes)) {
    stop(sprintf("Could not parse byte size for '%s'.", package_path), call. = FALSE)
  }

  description_path <- file.path(package_path, "DESCRIPTION")
  package_version <- read.dcf(description_path)[1, "Version"]

  list(
    path = package_path,
    size_human = size_human,
    size_bytes = size_bytes,
    version = package_version
  )
}

install_from_cran <- function(package_name, cran_mirror, library_dir) {
  utils::install.packages(
    package_name,
    repos = cran_mirror,
    lib = library_dir,
    quiet = TRUE
  )

  measure_installed_package(package_name, library_dir)
}

ensure_remotes <- function(cran_mirror) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes", repos = cran_mirror, quiet = TRUE)
  }
}

install_from_github <- function(package_name, github_repo, github_ref, library_dir, cran_mirror) {
  ensure_remotes(cran_mirror)

  remotes::install_github(
    repo = github_repo,
    ref = github_ref,
    lib = library_dir,
    upgrade = "never",
    quiet = TRUE,
    dependencies = TRUE
  )

  measure_installed_package(package_name, library_dir)
}

format_size_change <- function(size_bytes) {
  sign <- if (size_bytes > 0) "+" else ""
  sprintf("%s%.2f MB", sign, size_bytes / 1024 ^ 2)
}

format_percent_change <- function(size_bytes, baseline_bytes) {
  if (baseline_bytes == 0) {
    return("N/A")
  }

  sign <- if (size_bytes > 0) "+" else ""
  sprintf("%s%.1f%%", sign, (size_bytes / baseline_bytes) * 100)
}

write_reports <- function(output_dir, package_name, cran_result, github_result, github_repo, github_ref) {
  difference_bytes <- github_result$size_bytes - cran_result$size_bytes

  summary_lines <- c(
    sprintf("# %s package size comparison", package_name),
    "",
    sprintf("Generated: %s UTC", format(Sys.time(), tz = "UTC", usetz = FALSE)),
    "",
    "This report intentionally uses a very simple comparison:",
    "",
    sprintf("1. `install.packages(\"%s\")`", package_name),
    sprintf("2. `du -sh $(Rscript -e 'cat(find.package(\"%s\"))')`", package_name),
    sprintf("3. `remotes::install_github(\"%s\", ref = \"%s\")`", github_repo, github_ref),
    sprintf("4. `du -sh $(Rscript -e 'cat(find.package(\"%s\"))')`", package_name),
    "",
    "## Summary",
    "",
    "| Source | Version | Installed size | Path |",
    "| --- | --- | ---: | --- |",
    sprintf("| CRAN | %s | %s | `%s` |", cran_result$version, cran_result$size_human, cran_result$path),
    sprintf("| GitHub (`%s`) | %s | %s | `%s` |", github_ref, github_result$version, github_result$size_human, github_result$path),
    sprintf("| Difference (GitHub - CRAN) | - | %s | - |", format_size_change(difference_bytes)),
    sprintf("| Percent difference | - | %s | - |", format_percent_change(difference_bytes, cran_result$size_bytes))
  )

  report_lines <- c(
    summary_lines,
    "",
    "## Raw measurements",
    "",
    "| Source | `du -sh` | Bytes (`du -sb`) |",
    "| --- | ---: | ---: |",
    sprintf("| CRAN | %s | %s |", cran_result$size_human, format(cran_result$size_bytes, big.mark = ",", scientific = FALSE, trim = TRUE)),
    sprintf("| GitHub (`%s`) | %s | %s |", github_ref, github_result$size_human, format(github_result$size_bytes, big.mark = ",", scientific = FALSE, trim = TRUE)),
    "",
    "## Notes",
    "",
    "- Each install uses its own temporary library directory so the measurements stay isolated.",
    sprintf("- The GitHub install is taken from `%s` at ref `%s`.", github_repo, github_ref),
    "- The summary above is intended to be easy to read directly from the workflow run."
  )

  writeLines(summary_lines, con = file.path(output_dir, "package-size-summary.md"))
  writeLines(report_lines, con = file.path(output_dir, "package-size-report.md"))
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(arguments$output_dir, mustWork = FALSE)
cran_library <- tempfile("cran-library-")
github_library <- tempfile("github-library-")
dir.create(cran_library, recursive = TRUE, showWarnings = FALSE)
dir.create(github_library, recursive = TRUE, showWarnings = FALSE)

cran_result <- install_from_cran(
  package_name = arguments$package_name,
  cran_mirror = arguments$cran_mirror,
  library_dir = cran_library
)

github_result <- install_from_github(
  package_name = arguments$package_name,
  github_repo = arguments$github_repo,
  github_ref = arguments$github_ref,
  library_dir = github_library,
  cran_mirror = arguments$cran_mirror
)

write_reports(
  output_dir = output_dir,
  package_name = arguments$package_name,
  cran_result = cran_result,
  github_result = github_result,
  github_repo = arguments$github_repo,
  github_ref = arguments$github_ref
)
