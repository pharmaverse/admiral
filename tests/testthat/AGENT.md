# Admiral Unit Testing Guidelines for AI Assistants

Context for AI assistants when working with admiral unit tests in `tests/testthat/`.

**Auto-generated:** 2026-02-18 18:44:07
**Source:** admiraldev unit testing guidance vignette

---

# Admiral Unit Testing Guidelines

**Note:** Could not download latest admiraldev content. Using essential guidelines.
**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html

## Test Structure and Organization

### File Naming and Organization
- Match function names: `test-derive_var_example.R`
- One test file per main function (can include helpers)
- Use descriptive test names explaining the scenario

### Admiral Test Requirements

1. **Happy path** - basic functionality works
2. **Error conditions** - invalid inputs fail appropriately
3. **Edge cases** - empty data, boundary conditions
4. **Argument validation** - all parameters properly validated
5. **Output structure** - correct columns, types, ungrouped

### Test Data Best Practices

- Use `tribble()` for test data creation
- Test that output is ungrouped with `expect_false(is.grouped_df(result))`
- Test meaningful error messages
- Aim for 80%+ test coverage

---

*For complete guidance: https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html*
