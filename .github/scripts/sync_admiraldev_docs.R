#!/usr/bin/env Rscript

#' Sync admiraldev Documentation for GitHub Copilot
#'
#' This script copies key vignettes from the admiraldev package into
#' GitHub Copilot's instructions directory.
#'
#' Vignettes included:
#' - Programming Strategy
#' - Writing Vignettes
#' - Unit Test Guidance

# Suppress package startup messages
suppressPackageStartupMessages({
  library(fs)
  library(glue)
})

# Determine project root (works from any subdirectory)
if (require("here", quietly = TRUE)) {
  project_root <- here::here()
} else {
  # Fallback: find .git directory or use current directory
  project_root <- getwd()
  while (!dir_exists(file.path(project_root, ".git")) &&
    project_root != dirname(project_root)) {
    project_root <- dirname(project_root)
  }
}

# Configuration (relative to project root)
output_file <- file.path(project_root, ".github", "copilot-instructions.md")
log_dir <- file.path(project_root, "logs")
use_github <- TRUE # Set to TRUE to pull from GitHub, FALSE to use local package

cat(glue("Working from project root: {project_root}\n\n"))

# Create directories if they don't exist
dir_create(dirname(output_file), recurse = TRUE)
dir_create(log_dir, recurse = TRUE)

# Set up logging
log_file <- file.path(log_dir, paste0("sync_admiraldev_", Sys.Date(), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=================================\n")
cat("admiraldev Documentation Sync\n")
cat("Started:", format(Sys.time()), "\n")
cat("=================================\n\n")

# Main execution wrapped in tryCatch
tryCatch({
  # Define vignettes to copy with their display names
  vignettes_info <- list(
    list(
      file = "programming_strategy.Rmd",
      title = "Programming Strategy",
      description = "Core programming principles and strategies for admiral packages"
    ),
    list(
      file = "writing_vignettes.Rmd",
      title = "Writing Vignettes",
      description = "Guidelines for creating comprehensive and useful vignettes"
    ),
    list(
      file = "unit_test_guidance.Rmd",
      title = "Unit Test Guidance",
      description = "Best practices for writing unit tests in admiral packages"
    )
  )

  if (use_github) {
    cat("Using GitHub source for vignettes...\n\n")

    # GitHub raw URL base
    github_base <- "https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes"

    # Function to download vignette from GitHub
    download_vignette <- function(vignette_info) {
      file_name <- vignette_info$file
      title <- vignette_info$title
      description <- vignette_info$description

      url <- file.path(github_base, file_name)

      cat(glue("✓ Downloading: {file_name} from GitHub\n"))

      tryCatch(
        {
          # Download content
          content <- readLines(url, warn = FALSE)

          # Create formatted section
          formatted <- glue("

# {title}

**Description:** {description}
**Source File:** `{file_name}`
**Source:** GitHub (pharmaverse/admiraldev)
**URL:** {url}

---

{paste(content, collapse = '\n')}

---

")
          return(formatted)
        },
        error = function(e) {
          warning(glue("✗ Failed to download {file_name}: {e$message}"))
          return(NULL)
        }
      )
    }

    admiraldev_version <- "GitHub main branch"
    admiraldev_path <- "https://github.com/pharmaverse/admiraldev"
  } else {
    cat("Using local package installation...\n\n")

    # Try to find admiraldev package location
    admiraldev_path <- system.file(package = "admiraldev")

    if (admiraldev_path == "") {
      stop(
        "admiraldev package not found. Please install it with:\n",
        "  install.packages('admiraldev', dependencies = TRUE, build_vignettes = TRUE)\n",
        "or from GitHub:\n",
        "  remotes::install_github('pharmaverse/admiraldev', build_vignettes = TRUE)\n",
        "Or set use_github = TRUE in this script to download from GitHub directly."
      )
    }

    cat(glue("✓ Found admiraldev at: {admiraldev_path}\n"))

    # Check multiple possible vignette locations
    possible_dirs <- c(
      file.path(admiraldev_path, "doc"),
      file.path(admiraldev_path, "vignettes"),
      file.path(admiraldev_path, "inst", "doc")
    )

    vignettes_dir <- NULL
    for (dir in possible_dirs) {
      if (dir_exists(dir)) {
        # Check if it actually contains vignette files
        files <- list.files(dir, pattern = "\\.Rmd$|\\.html$")
        if (length(files) > 0) {
          vignettes_dir <- dir
          cat(glue("✓ Found vignettes in: {dir}\n"))
          cat(glue("  Files found: {paste(head(files, 3), collapse = ', ')}\n"))
          break
        }
      }
    }

    if (is.null(vignettes_dir)) {
      stop(glue(
        "Vignettes directory not found. Tried:\n",
        "  {paste(possible_dirs, collapse = '\n  ')}\n\n",
        "The package may have been installed without vignettes.\n",
        "Try reinstalling with:\n",
        "  install.packages('admiraldev', dependencies = TRUE, build_vignettes = TRUE)\n",
        "Or set use_github = TRUE in this script."
      ))
    }

    cat("\n")

    # Function to process and format vignette content from local files
    download_vignette <- function(vignette_info) {
      file_name <- vignette_info$file
      title <- vignette_info$title
      description <- vignette_info$description

      # Try .Rmd first, then .html as fallback
      rmd_path <- file.path(vignettes_dir, file_name)
      html_path <- file.path(vignettes_dir, sub("\\.Rmd$", ".html", file_name))

      if (file_exists(rmd_path)) {
        cat(glue("✓ Processing: {file_name}\n"))
        content <- readLines(rmd_path, warn = FALSE)
        source_path <- rmd_path
        source_type <- "R Markdown"
      } else if (file_exists(html_path)) {
        cat(glue("⚠ Using HTML version for: {file_name}\n"))
        content <- readLines(html_path, warn = FALSE)
        source_path <- html_path
        source_type <- "HTML"
      } else {
        warning(glue("✗ Vignette not found: {file_name}"))
        return(NULL)
      }

      # Create formatted section
      formatted <- glue("

# {title}

**Description:** {description}
**Source File:** `{basename(source_path)}`
**Format:** {source_type}

---

{paste(content, collapse = '\n')}

---

")

      return(formatted)
    }

    admiraldev_version <- as.character(packageVersion("admiraldev"))
  }

  # Create table of contents
  toc_items <- character(length(vignettes_info))
  for (i in seq_along(vignettes_info)) {
    title <- vignettes_info[[i]]$title
    anchor <- tolower(gsub("[^a-z0-9]+", "-", title))
    toc_items[i] <- paste0(i, ". [", title, "](#", anchor, ")")
  }
  toc_text <- paste(toc_items, collapse = "\n")

  # Create header for the combined file
  header <- glue("
# GitHub Copilot Instructions - admiraldev Guidelines

This file provides context to GitHub Copilot about admiraldev programming standards
and best practices. Copilot will automatically reference these guidelines when
providing code suggestions in this repository.

**Auto-generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}
**admiraldev Version:** {admiraldev_version}
**Source:** {admiraldev_path}

⚠️ **DO NOT EDIT MANUALLY** - Run `source('scripts/sync_admiraldev_docs.R')` to update

---

## Purpose

These guidelines ensure that code in this repository follows admiral ecosystem standards:

- ✅ Consistent programming patterns across admiral packages
- ✅ Proper function design and documentation
- ✅ High-quality vignettes that help users
- ✅ Comprehensive unit tests with good coverage
- ✅ Code that is maintainable and readable

GitHub Copilot will use these guidelines to provide better suggestions that align
with admiral best practices.

---

## Table of Contents

{toc_text}

---

")

  # Combine all content
  cat("Processing vignettes...\n")
  combined_content <- header

  processed_count <- 0
  for (vignette_info in vignettes_info) {
    vignette_content <- download_vignette(vignette_info)
    if (!is.null(vignette_content)) {
      combined_content <- paste0(combined_content, vignette_content)
      processed_count <- processed_count + 1
    }
  }

  # Add footer
  footer <- glue("

---

## Keeping This Document Updated

**Generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}
**admiraldev Version:** {admiraldev_version}
**Vignettes Included:** {processed_count} of {length(vignettes_info)}
**Output File:** `{output_file}`

### Manual Update
```r
source('scripts/sync_admiraldev_docs.R')
```

### Automatic Updates

Consider setting up a monthly cron job or GitHub Action to keep these
guidelines synchronized with the latest admiraldev version.

See `scripts/sync_admiraldev_docs.R` for setup instructions.

---

## How GitHub Copilot Uses This File

GitHub Copilot automatically reads files in `.github/` directories and uses them
as context when providing code suggestions. By keeping admiraldev guidelines here:

1. **Copilot suggests code** that follows admiral conventions
2. **Function names and patterns** match admiral ecosystem standards
3. **Documentation style** aligns with admiral expectations
4. **Test structures** follow unit test guidance
5. **Vignette examples** use recommended patterns

You don't need to do anything special - just having this file in `.github/`
is enough for Copilot to use it.

---

*This is an automated sync of admiraldev documentation for GitHub Copilot context.*
")

  combined_content <- paste0(combined_content, footer)

  # Write to output file
  writeLines(combined_content, output_file)

  # Summary
  cat("\n=================================\n")
  cat("✓ Sync completed successfully!\n")
  cat("=================================\n\n")
  cat(glue("📄 Output file: {output_file}\n"))
  cat(glue("📊 Processed: {processed_count}/{length(vignettes_info)} vignettes\n"))
  cat(glue("📝 Log file: {log_file}\n"))
  cat(glue("⏰ Completed: {format(Sys.time())}\n\n"))

  cat("GitHub Copilot will now use these guidelines automatically!\n")
  cat("No need to open the file - Copilot reads .github/ automatically.\n\n")

  # Return success status
  invisible(TRUE)
}, error = function(e) {
  cat("\n=================================\n")
  cat("✗ ERROR during sync\n")
  cat("=================================\n\n")
  cat("Error message:\n")
  cat(conditionMessage(e), "\n\n")
  cat("Traceback:\n")
  print(sys.calls())

  # Return failure status
  invisible(FALSE)
}, finally = {
  # Clean up logging
  sink(type = "message")
  sink(type = "output")
  close(log_con)
})
