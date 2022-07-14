# GitHub Actions Workflows

This directory contains GitHub Actions Continuous integration/continuous delivery (CI/CD) workflows, most of which are also used by `admiral`'s extensions.
Workflows defined here are responsible for assuring high package quality standards without compromising performance, security, or reproducibility.

## What these workflows do?

Most workflows have a `BEGIN boilderplate steps` and `END boilderplate steps` section within them which define some standard steps required for installing system dependencies, R version and R packages which serve as dependencies for the package.

The underlying mechanisms for installing R and Pandoc are defined in [`r-lib/actions`][r-lib-actions], while the installation of system dependencies and R package dependencies is managed via the [Staged Dependencies GitHub Action][sd-action]. The latter is used in conjunction with the [`staged_dependencies.yaml`](../../staged_dependencies.yaml) file in order to install dependencies that are in the _same stage of development_ as the current package. You can read more about how it works [here](sd-repo). Note that the latter is not necesary for this workflow to work and is completely optional.

Following the installation of system dependencies, R, and package dependencies, each workflow does something different.

### [`check-templates.yml`](check-templates.yml)

This workflow checks for issues within template scripts.  For example, in admiral package there are several template scripts with admiral-based functions showing how to build certain ADaM datasets.   As we update the admiral functions, we want to make sure these template scripts execute appropriately.  Functions in the template scripts that are deprecated or used inappropriately will cause this workflow to fail. 

### [`code-coverage.yml`](code-coverage.yml)

This workflow measures code coverage for unit tests and reports the code coverage as a percentage of the _total number of lines covered by unit tests_ vs. the _total number of lines in the codebase_.

The [`covr`][covr] R package is used to calculate the coverage.

Report summaries and badges for coverage are generated using a series of other GitHub Actions.

### [`links.yml`](links.yml)

This workflow checks whether URLs embedded in code and documentation are valid. Invalid URLs reults in workflow failures. This workflow uses [`lychee`][lychee] to detect broken links.   Occasionally this check will detect false positives of urls that look like urls.  To remedy, please add this false positive to the `.lycheeignore` file.   

### [`lintr.yml`](lintr.yml)

Static code analysis is performed by this workflow, which in turn uses the [`lintr`][lintr] R package.

Any [`.lintr`](../../.lintr) configurations in the repository will be by this workflow.

### [`man-pages.yml`](man-pages.yml)

This workflow checks if the manual pages in the `man/` directory of the package are up-to-date with ROxygen comments in the code.

Workflow failures indicate that the manual pages are not up-to-date with ROxygen comments, and corrective actions are provided in the workflow log.

### [`pkgdown.yml`](pkgdown.yml)

Documentation for the R package is generated via this workflow. This workflow uses the [`pkgdown`](pkgdown) framework to generate documentation in HTML, and the HTML pages are deployed to the `gh-pages` branch.

Moreover, an additional `Versions` dropdown is generated via the [][multi-version-docs] GitHub Action, so that an end user can view multiple versions of the documentation for the package.

### [`r-cmd-check.yml`](r-cmd-check.yml)

This workflow performs `R CND check` for the package. Failed workflows are typically indicative of problems encountered during the check, and therefore an indication that the package does not meet quality standards.

### [`r-pkg-validation.yml`](r-pkg-validation.yml)

When a new release of the package is made, this workflow executes to create a validation report via [theValidatoR][validation]. The PDF report is then attached to the release within GitHub.

### [`readme-render.yml`](readme-render.yml)

If your codebase uses a [`README.Rmd` file](../../README.Rmd), then this workflow will automatically render a `README.md` and commit it to your branch.

### [`spellcheck.yml`](spellcheck.yml)

Spellchecks are performed by this workflow, and the [`spelling`][spelling] R package is used to detect spelling mistakes. Failed workflows typically indicate misspelled words.  In the `inst/WORDLIST` file, you can add words and or acronyms that you want the spell check to ignore, for example occds is not an English word but a common acronym used within Pharma.  The workflow will flag this until a user adds it to the `inst/WORDLIST`.

### [`style.yml`](style.yml`)

Code style is enforced via the [`styler`][styler] R package. Custom style configurations, if any, will be honored by this workflow.

Failed workflows are indicative of unstyled code.

## How to use these workflows?

### Reuse (recommended)

You could add just _one_ file called `.github/workflows/common.yml` to directly import these workflows while receiving the latest updates and enhancements, given that the workflows defined in this repository enable reusability via the [`workflow_call`][workflow_call] GitHub Actions event.

The contents of the `.github/workflows/common.yml` file are available in the [`common.yml.inactive`](common.yml.inactive) file in this repository.

### Copy as-is (not recommended)

Alternatively, if you want a high level of customization, you could simply copy the workflows as-is from this repository to your repository.

<!-- Begin links -->
[r-lib-actions]: https://github.com/r-lib/actions
[sd-action]: https://github.com/marketplace/actions/staged-dependencies-action
[sd-repo]: https://github.com/openpharma/staged.dependencies
[lychee]: https://github.com/lycheeverse/lychee
[covr]: https://covr.r-lib.org/
[lintr]: https://lintr.r-lib.org/
[pkgdown]: https://pkgdown.r-lib.org/
[multi-version-docs]: https://github.com/marketplace/actions/r-pkgdown-multi-version-docs
[validation]: https://github.com/marketplace/actions/r-package-validation-report
[spelling]: https://docs.ropensci.org/spelling/
[styler]: https://styler.r-lib.org/
[workflow_call]: https://docs.github.com/en/actions/using-workflows/reusing-workflows
<!-- End links -->
