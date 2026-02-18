# Admiral Development Guidelines for AI Assistants

This file provides context for AI coding assistants (GitHub Copilot, Gemini, Claude, etc.) about admiral ecosystem standards and best practices.

**Auto-generated:** 2026-02-18 18:44:07
**Source:** admiraldev package vignettes
**Update script:** `source('.github/scripts/sync_admiraldev_copilot.R')`

## Purpose

Help AI assistants provide admiral-compliant code suggestions:
- ✅ Consistent function naming and patterns
- ✅ Proper argument handling and validation
- ✅ CDISC/ADaM context and conventions
- ✅ Testing and documentation standards

---

# Admiral Programming Strategy

**Note:** Could not download latest admiraldev content. Using essential guidelines.
**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html

## Function Design Principles

- **Modularity:** All code follows a modular approach with clearly separated steps
- **Avoid Copy and Paste:** Similar code should be put into separate functions
- **Checks:** Meaningful error messages with clear reference to failing input
- **Flexibility:** Functions should be as flexible as possible without reducing usability

## Function Naming Convention

Function names should start with a verb and use snake_case:

- `derive_*` - Add rows/columns to datasets
- `derive_var_*` - Add single variable
- `derive_vars_*` - Add multiple variables
- `compute_*` - Vector operations, return vectors
- `assert_*` / `warn_*` - Input validation
- `filter_*` - Filter observations

---

## Testing Guidelines

For admiral unit testing guidelines, see `tests/testthat/AGENT.md`

## Full Documentation

- **Programming Strategy:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html
- **Unit Test Guidance:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html
- **Admiral Website:** https://pharmaverse.github.io/admiral/

---

*This file helps AI assistants understand admiral patterns for better code suggestions.*
