# admiral guidelines loaded

parse_args <- function(args) {
  defaults <- list(
    repo = getwd(),
    output_dir = file.path(getwd(), "reports"),
    label = "main",
    cran_mirror = "https://cran.r-project.org"
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

run_command <- function(command, args, working_directory = NULL) {
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)

  if (!is.null(working_directory)) {
    setwd(working_directory)
  }

  output <- system2(
    command = command,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste(
        c(
          sprintf("Command failed: %s %s", command, paste(args, collapse = " ")),
          output
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  invisible(output)
}

build_source_tarball <- function(repo_dir, build_dir) {
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

  package_dir <- normalizePath(repo_dir, winslash = "/", mustWork = TRUE)
  package_name <- basename(package_dir)
  parent_dir <- dirname(package_dir)

  run_command(
    command = file.path(R.home("bin"), "R"),
    args = c("CMD", "build", package_name),
    working_directory = parent_dir
  )

  tarballs <- list.files(parent_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  package_tarballs <- tarballs[basename(tarballs) %in% sprintf("%s_*.tar.gz", package_name)]

  if (length(package_tarballs) == 0) {
    package_tarballs <- tarballs[grepl(sprintf("^%s_.*\\.tar\\.gz$", package_name), basename(tarballs))]
  }

  if (length(package_tarballs) == 0) {
    stop("No source tarball was created by `R CMD build`.", call. = FALSE)
  }

  latest_tarball <- package_tarballs[[which.max(file.info(package_tarballs)$mtime)]]
  destination <- file.path(build_dir, basename(latest_tarball))
  success <- file.copy(latest_tarball, destination, overwrite = TRUE)

  if (!isTRUE(success)) {
    stop("Failed to copy the built source tarball to the build directory.", call. = FALSE)
  }

  destination
}

download_cran_tarball <- function(package_name, cran_mirror, download_dir) {
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

  cran_db <- utils::available.packages(
    repos = cran_mirror,
    type = "source"
  )

  if (!package_name %in% rownames(cran_db)) {
    stop(
      sprintf("Package '%s' was not found in the CRAN package database.", package_name),
      call. = FALSE
    )
  }

  cran_version <- cran_db[package_name, "Version"]
  cran_url <- sprintf(
    "%s/%s_%s.tar.gz",
    utils::contrib.url(cran_mirror, type = "source"),
    package_name,
    cran_version
  )
  destination <- file.path(download_dir, basename(cran_url))

  utils::download.file(
    url = cran_url,
    destfile = destination,
    mode = "wb",
    quiet = TRUE
  )

  list(path = destination, version = cran_version)
}

collect_tarball_inventory <- function(tarball) {
  extract_dir <- tempfile("package-size-")
  dir.create(extract_dir)
  utils::untar(tarball, exdir = extract_dir)

  package_roots <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
  if (length(package_roots) == 0) {
    stop(
      sprintf("No package contents were extracted from '%s'.", tarball),
      call. = FALSE
    )
  }

  package_root <- package_roots[[1]]
  files <- list.files(
    package_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  normalized_root <- normalizePath(package_root, winslash = "/", mustWork = TRUE)
  normalized_files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  relative_paths <- sub(paste0("^", normalized_root, "/"), "", normalized_files)
  file_sizes <- unname(file.info(files)$size)

  inventory <- data.frame(
    path = relative_paths,
    size_bytes = file_sizes,
    size_mb = round(file_sizes / 1024 ^ 2, 6),
    stringsAsFactors = FALSE
  )

  inventory <- inventory[order(-inventory$size_bytes, inventory$path), ]
  inventory$rank <- seq_len(nrow(inventory))
  inventory[, c("rank", "path", "size_bytes", "size_mb")]
}

format_mb <- function(size_bytes) {
  sprintf("%.2f", size_bytes / 1024 ^ 2)
}

write_markdown_report <- function(
    output_file,
    package_name,
    cran_version,
    cran_tarball_size_bytes,
    cran_extracted_size_bytes,
    development_version,
    development_label,
    development_tarball_size_bytes,
    development_extracted_size_bytes,
    development_inventory
) {
  tarball_size_difference_bytes <- development_tarball_size_bytes - cran_tarball_size_bytes
  tarball_percent_difference <- if (cran_tarball_size_bytes == 0) {
    NA_real_
  } else {
    (tarball_size_difference_bytes / cran_tarball_size_bytes) * 100
  }

  extracted_size_difference_bytes <- development_extracted_size_bytes - cran_extracted_size_bytes
  extracted_percent_difference <- if (cran_extracted_size_bytes == 0) {
    NA_real_
  } else {
    (extracted_size_difference_bytes / cran_extracted_size_bytes) * 100
  }

  top_files <- utils::head(development_inventory, n = 10)
  top_file_lines <- paste0(
    "| ",
    top_files$rank,
    " | `",
    top_files$path,
    "` | ",
    format(top_files$size_bytes, big.mark = ",", scientific = FALSE, trim = TRUE),
    " | ",
    sprintf("%.2f", top_files$size_mb),
    " |"
  )

  report_lines <- c(
    sprintf("# %s package size comparison", package_name),
    "",
    sprintf("Generated: %s UTC", format(Sys.time(), tz = "UTC", usetz = FALSE)),
    "",
    "## Source tarball size (CRAN-relevant)",
    "",
    "| Source | Version | Size (MB) |",
    "| --- | --- | ---: |",
    sprintf("| CRAN | %s | %s |", cran_version, format_mb(cran_tarball_size_bytes)),
    sprintf(
      "| Development (%s) | %s | %s |",
      development_label,
      development_version,
      format_mb(development_tarball_size_bytes)
    ),
    sprintf("| Difference (development - CRAN) | - | %s |", format_mb(tarball_size_difference_bytes)),
    sprintf(
      "| Percent difference | - | %s%% |",
      if (is.na(tarball_percent_difference)) "N/A" else sprintf("%.2f", tarball_percent_difference)
    ),
    "",
    "## Extracted package contents size (diagnostic)",
    "",
    "| Source | Version | Size (MB) |",
    "| --- | --- | ---: |",
    sprintf("| CRAN | %s | %s |", cran_version, format_mb(cran_extracted_size_bytes)),
    sprintf(
      "| Development (%s) | %s | %s |",
      development_label,
      development_version,
      format_mb(development_extracted_size_bytes)
    ),
    sprintf("| Difference (development - CRAN) | - | %s |", format_mb(extracted_size_difference_bytes)),
    sprintf(
      "| Percent difference | - | %s%% |",
      if (is.na(extracted_percent_difference)) "N/A" else sprintf("%.2f", extracted_percent_difference)
    ),
    "",
    sprintf(
      "The development tarball contains %s files after `R CMD build` packaging.",
      format(nrow(development_inventory), big.mark = ",", trim = TRUE)
    ),
    "",
    "## Largest files in the development package",
    "",
    "| Rank | File | Size (bytes) | Size (MB) |",
    "| ---: | --- | ---: | ---: |",
    top_file_lines
  )

  writeLines(report_lines, con = output_file)
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
repo_dir <- normalizePath(arguments$repo, mustWork = TRUE)
dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(arguments$output_dir, mustWork = FALSE)

description <- read.dcf(file.path(repo_dir, "DESCRIPTION"))
package_name <- description[1, "Package"]
development_version <- description[1, "Version"]

development_tarball <- build_source_tarball(
  repo_dir = repo_dir,
  build_dir = file.path(tempdir(), "development-build")
)
cran_download <- download_cran_tarball(
  package_name = package_name,
  cran_mirror = arguments$cran_mirror,
  download_dir = file.path(tempdir(), "cran-download")
)

development_inventory <- collect_tarball_inventory(development_tarball)
cran_inventory <- collect_tarball_inventory(cran_download$path)

development_tarball_size_bytes <- unname(file.info(development_tarball)$size)
cran_tarball_size_bytes <- unname(file.info(cran_download$path)$size)
development_extracted_size_bytes <- sum(development_inventory$size_bytes)
cran_extracted_size_bytes <- sum(cran_inventory$size_bytes)

csv_output <- file.path(output_dir, "package-size-development-files.csv")
markdown_output <- file.path(output_dir, "package-size-report.md")

utils::write.csv(development_inventory, csv_output, row.names = FALSE)

write_markdown_report(
  output_file = markdown_output,
  package_name = package_name,
  cran_version = cran_download$version,
  cran_tarball_size_bytes = cran_tarball_size_bytes,
  cran_extracted_size_bytes = cran_extracted_size_bytes,
  development_version = development_version,
  development_label = arguments$label,
  development_tarball_size_bytes = development_tarball_size_bytes,
  development_extracted_size_bytes = development_extracted_size_bytes,
  development_inventory = development_inventory
)
