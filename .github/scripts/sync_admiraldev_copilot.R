#!/usr/bin/env Rscript

#' Sync admiraldev Documentation for GitHub Copilot (v1.0)
#'
#' Minimal implementation for admiral 1.4.0: downloads key admiraldev vignettes
#' and creates Copilot instructions with essential admiral programming context.
#'
#' Handles SSL/network issues gracefully with fallback content.

# Load required packages
suppressPackageStartupMessages({
  library(glue)
})

# Configuration
project_root <- if (file.exists(".git")) "." else ".."
output_file <- file.path(project_root, ".github", "copilot-instructions.md")
github_raw_base <- "https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes"

# Create .github directory if needed
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

cat("Syncing admiraldev guidelines for GitHub Copilot...\n")

# Define the 2 most important vignettes for admiral development
vignettes <- list(
  list(
    file = "programming_strategy.Rmd",
    title = "Programming Strategy",
    description = "Core programming principles and strategies for admiral packages"
  ),
  list(
    file = "unit_test_guidance.Rmd",
    title = "Unit Test Guidance",
    description = "Best practices for writing unit tests in admiral packages"
  )
)

# Improved download function with better error handling
process_vignette <- function(vignette_info) {
  url <- file.path(github_raw_base, vignette_info$file)

  cat(glue("• {vignette_info$title}..."))

  # Try multiple download methods
  content <- NULL

  # Method 1: readLines (simplest)
  tryCatch(
    {
      content <- readLines(url, warn = FALSE)
      cat(" ✓ (readLines)\n")
    },
    error = function(e) {
      # Method 2: download.file + readLines (sometimes works better with SSL)
      tryCatch(
        {
          temp_file <- tempfile(fileext = ".Rmd")
          download.file(url, temp_file, quiet = TRUE, method = "auto")
          content <<- readLines(temp_file, warn = FALSE)
          unlink(temp_file)
          cat(" ✓ (download.file)\n")
        },
        error = function(e2) {
          # Method 3: curl if available
          if (Sys.which("curl") != "") {
            tryCatch(
              {
                temp_file <- tempfile(fileext = ".Rmd")
                system2("curl", c("-s", "-o", temp_file, url), stdout = FALSE, stderr = FALSE)
                if (file.exists(temp_file) && file.size(temp_file) > 0) {
                  content <<- readLines(temp_file, warn = FALSE)
                  unlink(temp_file)
                  cat(" ✓ (curl)\n")
                } else {
                  stop("curl failed")
                }
              },
              error = function(e3) {
                cat(" ✗ (SSL/network error)\n")
                warning(glue("Failed to download {vignette_info$file}: {e$message}"), call. = FALSE)
              }
            )
          } else {
            cat(" ✗ (SSL/network error)\n")
            warning(glue("Failed to download {vignette_info$file}: {e$message}"), call. = FALSE)
          }
        }
      )
    }
  )

  if (!is.null(content) && length(content) > 0) {
    return(glue("
# {vignette_info$title}

**Description:** {vignette_info$description}
**Source:** admiraldev vignette `{vignette_info$file}`
**URL:** {url}

---

{paste(content, collapse = '\n')}

---
"))
  } else {
    # Return fallback content if download failed
    return(create_fallback_content(vignette_info))
  }
}

# Fallback content for when downloads fail
create_fallback_content <- function(vignette_info) {
  if (vignette_info$file == "programming_strategy.Rmd") {
    return(glue("
# {vignette_info$title}

**Description:** {vignette_info$description}
**Source:** admiraldev vignette `{vignette_info$file}` (fallback content)
**Note:** Could not download from GitHub. Using essential guidelines.

---


---
"))
  } else if (vignette_info$file == "unit_test_guidance.Rmd") {
    return(glue("
# {vignette_info$title}

**Description:** {vignette_info$description}
**Source:** admiraldev vignette `{vignette_info$file}` (fallback content)
**Note:** Could not download from GitHub. Using essential guidelines.

---

**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html

---
"))
  } else {
    return(glue("
# {vignette_info$title}

**Description:** {vignette_info$description}
**Source:** admiraldev vignette `{vignette_info$file}` (fallback - could not download)

Please see: https://pharmaverse.github.io/admiraldev/articles/{gsub('\\.Rmd$', '.html', vignette_info$file)}

---
"))
  }
}

# Create header
header <- glue("# GitHub Copilot Instructions - admiral Development

**Auto-generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}
**Source:** admiraldev package vignettes

This file provides GitHub Copilot with context about admiral programming
standards and best practices. Copilot automatically uses these guidelines
when providing code suggestions in this repository.

⚠️ **DO NOT EDIT MANUALLY** - Run `source('.github/scripts/sync_admiraldev_copilot.R')` to update

## Purpose

Ensure code follows admiral ecosystem standards:
- Consistent programming patterns across admiral functions
- Proper function design and documentation
- Comprehensive unit tests with good coverage
- Maintainable and readable code

GitHub Copilot will use these guidelines to provide better suggestions that align
with admiral best practices.

---
")

# Process all vignettes
cat("Downloading vignettes from admiraldev:\n")
combined_content <- header

processed_count <- 0
for (vignette in vignettes) {
  vignette_content <- process_vignette(vignette)
  if (!is.null(vignette_content)) {
    combined_content <- paste0(combined_content, vignette_content)
    processed_count <- processed_count + 1
  }
}

# Add footer
footer <- glue("
---

## How to Update

```r
source('.github/scripts/sync_admiraldev_copilot.R')
```

## Network Issues?

If you see SSL or connection errors:
1. Try running from a different network
2. Use corporate VPN if available
3. The script includes fallback content for essential guidelines
4. Manual alternative: Copy content from https://pharmaverse.github.io/admiraldev/

## GitHub Copilot Integration

GitHub Copilot automatically reads files in `.github/` directories and uses them
as context when providing code suggestions. By keeping admiraldev guidelines here:

1. **Copilot suggests code** that follows admiral conventions
2. **Function names and patterns** match admiral ecosystem standards
3. **Documentation style** aligns with admiral expectations
4. **Test structures** follow unit test guidance

You don't need to do anything special - just having this file in `.github/`
is enough for Copilot to use it.

---

**Content Status:** Successfully processed {processed_count}/{length(vignettes)} vignettes
*Generated from admiraldev vignettes to provide GitHub Copilot with admiral context.*
")

combined_content <- paste0(combined_content, footer)

# Write output
writeLines(combined_content, output_file)

cat(glue("\n✓ Created: {output_file}\n"))
cat(glue("📊 Processed: {processed_count}/{length(vignettes)} vignettes\n"))

if (processed_count == length(vignettes)) {
  cat("🎉 All vignettes downloaded successfully!\n")
} else {
  cat("⚠️  Some downloads failed - using fallback content\n")
  cat("💡 Try running from a different network or with VPN\n")
}

cat("\nGitHub Copilot will now use admiral guidelines automatically!\n")
