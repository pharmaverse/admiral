# GitHub Copilot Instructions - admiral Development

**Auto-generated:** 2026-02-18 18:57:31
**Optimized for:** GitHub Copilot code completion
**Complete guidelines:** See `AGENT.md` files for full context

⚠️ **DO NOT EDIT MANUALLY** - Run `source('.github/scripts/sync_admiraldev_copilot.R')` to update

## Admiral Code Completion Patterns

### Function Names (verb_object_detail)
- `derive_var_base()` - Single variable derivation
- `derive_vars_merged()` - Multiple variables from merge
- `derive_param_tte()` - Parameter derivation
- `compute_age_years()` - Vector computation
- `assert_data_frame()` - Input validation

### Argument Patterns
Standard admiral function signature:
- `dataset` - Always first argument
- `by_vars = exprs(...)` - Grouping variables
- `order = exprs(...)` - Sorting expressions
- `new_var = VAR` - New variable (symbol)
- `filter = PARAM == "VALUE"` - Filtering expression

### Common Expressions
- `subject_keys = exprs(STUDYID, USUBJID)`
- `by_vars = exprs(USUBJID, PARAMCD)`
- `order = exprs(AVISITN, desc(ADY))`
- `filter = PARAMCD == "TEMP"`

### Test Patterns
- Use `tribble()` for test data
- `expect_true("DERIVED" %in% names(result))`
- `expect_false(is.grouped_df(result))`
- `expect_error(..., "Required variable.*missing")`

---

*This file is optimized for GitHub Copilot code completion. For complete admiral development guidelines, see `AGENT.md` and `tests/testthat/AGENT.md`.*
