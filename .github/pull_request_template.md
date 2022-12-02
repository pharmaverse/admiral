Thank you for your Pull Request! We have developed this task checklist from the [Development Process Guide](https://pharmaverse.github.io/admiraldev/devel/articles/development_process.html) to help with the final steps of the process. Completing the below tasks helps to ensure our reviewers can maximize their time on your code as well as making sure the admiral codebase remains robust and consistent.   

Please check off each taskbox as an acknowledgment that you completed the task or check off that it is not relevant to your Pull Request. This checklist is part of the Github Action workflows and the Pull Request will not be merged into the `devel` branch until you have checked off each task.

- [ ] Place Closes #<insert_issue_number> into the beginning of your Pull Request Title (Use Edit button in top-right if you need to update)
- [ ] Code is formatted according to the [tidyverse style guide](https://style.tidyverse.org/). Run `styler::style_file()` to style R and Rmd files
- [ ] Updated relevant unit tests or have written new unit tests, which should consider realistic data scenarios and edge cases, e.g. empty datasets, errors, boundary cases etc. - See [Unit Test Guide](https://pharmaverse.github.io/admiraldev/devel/articles/unit_test_guidance.html#tests-should-be-robust-to-cover-realistic-data-scenarios)
- [ ] If you removed/replaced any function and/or function parameters, did you fully follow the [deprecation guidance](https://pharmaverse.github.io/admiraldev/devel/articles/programming_strategy.html#deprecation)?
- [ ] Update to all relevant roxygen headers and examples, including keywords and families. Refer to the [categorization of functions](https://pharmaverse.github.io/admiraldev/devel/articles/programming_strategy.html#categorization-of-functions) to tag appropriate keyword/family.
- [ ] Run `devtools::document()` so all `.Rd` files in the `man` folder and the `NAMESPACE` file in the project root are updated appropriately
- [ ] Address any updates needed for vignettes and/or templates
- [ ] Update `NEWS.md` if the changes pertain to a user-facing function (i.e. it has an `@export` tag) or documentation aimed at users (rather than developers)
- [ ] Build admiral site `pkgdown::build_site()` and check that all affected examples are displayed correctly and that all new functions occur on the "[Reference](https://pharmaverse.github.io/admiral/reference/index.html)" page. 
- [ ] Address or fix all lintr warnings and errors - `lintr::lint_package()`
- [ ] Run `R CMD check` locally and address all errors and warnings - `devtools::check()`
- [ ] Link the issue in the Development Section on the right hand side.
- [ ] Address all merge conflicts and resolve appropriately 
- [ ] Pat yourself on the back for a job well done! Much love to your accomplishment!
