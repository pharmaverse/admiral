# admiral <img src="man/figures/logo.png" align="right" width="200" style="margin-left:50px;"/>

<!-- start badges -->
[![pharmaverse admiral Badge](http://pharmaverse.org/shields/admiral.svg)](https://pharmaverse.org)
[![CRAN status](https://www.r-pkg.org/badges/version/admiral)](https://CRAN.R-project.org/package=admiral)
![Test Coverage](https://raw.githubusercontent.com/pharmaverse/admiral/badges/main/test-coverage.svg)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/admiral?color=green)](https://cran.r-project.org/package=admiral)
[![admiral CI/CD Workflows](https://github.com/pharmaverse/admiral/actions/workflows/common.yml/badge.svg)](https://github.com/pharmaverse/admiral/actions/workflows/common.yml)

<!-- end badges -->
ADaM in R Asset Library

*Explore all the other packages in the [{admiral} ecosystem](https://pharmaverse.org/e2eclinical/adam/) to learn more about ADaM programming in R.*

## Purpose

To provide an open source, modularized toolbox that enables the pharmaceutical programming community
to develop ADaM datasets in R.

## Installation

The package is available from CRAN and can be installed with:

```r
install.packages("admiral")
```

To install the development version of the package from GitHub run:

```r
pak::pkg_install("pharmaverse/admiral", dependencies = TRUE)
```

## Cheat Sheet

<a href="https://github.com/pharmaverse/admiral/blob/main/inst/cheatsheet/admiral_cheatsheet.pdf"><img src="https://raw.githubusercontent.com/pharmaverse/admiral/main/inst/cheatsheet/cheatsheet_image.png" width="630" height="252" alt="Admiral Cheat Sheet Image"/></a>

## Release Schedule

The `{admiral}` family has several downstream and upstream dependencies and so releases are done in two Phases:

* Phase 1 release is for [{admiraldev}](https://pharmaverse.github.io/admiraldev/), [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/), and [{admiral}](https://pharmaverse.github.io/admiral/cran-release/) core
* Phase 2 release is extension packages, e.g. [{admiralonco}](https://pharmaverse.github.io/admiralonco/), [{admiralophtha}](https://pharmaverse.github.io/admiralophtha/), [{admiralvaccine}](https://pharmaverse.github.io/admiralvaccine/) and [{pharmaverseadam}](https://pharmaverse.github.io/pharmaverseadam/).

__NB:__ We strive for a regular 6 month release schedule for [{admiraldev}](https://pharmaverse.github.io/admiraldev/), [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/), and [{admiral}](https://pharmaverse.github.io/admiral/). Extension packages releases are on a content-basis and as such may be more infrequent than the below schedule shows, or there may even be ad-hoc releases.

| Release Schedule | Phase 1- Date and Packages                                               | Phase 2- Date and Packages                                             |
| ---------------- | ------------------------------------------------------------------------ | ---------------------------------------------------------------------- |
|                  |                                                                          |                                                                        |
| Q2 2026          | Mid-June 2026                                                            | End of June 2026                                                      | 
|                  | [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/)      | [{admiralonco}](https://pharmaverse.github.io/admiralonco/)            |
|                  | [{admiraldev}](https://pharmaverse.github.io/admiraldev/)                | [{admiralophtha}](https://pharmaverse.github.io/admiralophtha/)        |
|                  | [{admiral}](https://pharmaverse.github.io/admiral/)                      | [{admiralvaccine}](https://pharmaverse.github.io/admiralvaccine/)      |
|                  |                                                                          | [{admiralpeds}](https://pharmaverse.github.io/admiralpeds/)            |
|                  |                                                                          | [{admiralmetabolic}](https://pharmaverse.github.io/admiralmetabolic/)  |
|                  |                                                                          | [{admiralneuro}](https://pharmaverse.github.io/admiralneuro/)          |
|                  |                                                                          | [{pharmaverseadam}](https://pharmaverse.github.io/pharmaverseadam/)    |
|                  |                                                                          |                                                                        | 
| Q4 2026/Q1 2027  | Late December 2026/Early January 2027                                    |  Mid-January 2027                                                                      |
|                  | [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/)      | [{admiralonco}](https://pharmaverse.github.io/admiralonco/)            |
|                  | [{admiraldev}](https://pharmaverse.github.io/admiraldev/)                | [{admiralophtha}](https://pharmaverse.github.io/admiralophtha/)        |
|                  | [{admiral}](https://pharmaverse.github.io/admiral/)                      | [{admiralvaccine}](https://pharmaverse.github.io/admiralvaccine/)      |
|                  |                                                                          | [{admiralpeds}](https://pharmaverse.github.io/admiralpeds/)            |
|                  |                                                                          | [{admiralmetabolic}](https://pharmaverse.github.io/admiralmetabolic/)  |
|                  |                                                                          | [{admiralneuro}](https://pharmaverse.github.io/admiralneuro/)          |
|                  |                                                                          | [{pharmaverseadam}](https://pharmaverse.github.io/pharmaverseadam/)    |                                                                       

## Main Goal

Provide users with an open source, modularized toolbox with which to create ADaM datasets
in R. _As opposed to a "run one line and an ADaM appears" black-box solution or an attempt to
automate ADaM._

One of the key aspects of `{admiral}` is its development by the users for the users.
It gives an entry point for all to collaborate, co-create and contribute to a
harmonized approach of developing ADaMs in R across the pharmaceutical industry.

## Scope

To set expectations: It is not our target that `{admiral}` will ever provide all possible solutions
for all ADaM datasets outside of study specific needs. It depends on the user's collaboration
and contribution to help grow over time to an asset library that is robust, easy to use and
has an across-industry focus. We do not see a coverage of 100% of all ADaM derivations as ever
achievable---ADaM is endless.

We will provide:

* A toolbox of re-usable functions and utilities to create ADaM datasets using R scripts in a
  modular manner (an "opinionated" design strategy).
* Pharmaceutical communities and companies are encouraged to contribute to `{admiral}` following
  the provided programming strategy and modular approach
* Functions that are comprehensively documented and tested, including example calls---these are
  all listed in the [Reference section](https://pharmaverse.github.io/admiral/cran-release/reference/index.html).
* Vignettes on how to create ADSL, BDS and OCCDS datasets, including example scripts.
* Vignettes for ADaM dataset specific functionality (i.e. dictionary coding, date imputation, SMQs ...).

## The {admiral} Family of Packages

There are three types of packages in the `{admiral}` family:

* Core package---one package containing all core functions required to create ADaMs, usable by any company (i.e. general derivations, utility functions and checks for ADSL, OCCDS and BDS).
* TA (Therapeutic Area) package extensions---one package per TA with functions that are
  specific to algorithms and requirements for that particular TA (e.g. [`{admiralonco}`](https://pharmaverse.github.io/admiralonco/)).
* Company package extensions---specific needs and plug-ins for the company, such as access to metadata
  (e.g. `{admiralroche}` or `{admiralgsk}`).

## Related Packages

Related data packages include:

* [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/)---this contains test SDTM data sourced from the [CDISC pilot project](https://github.com/cdisc-org/sdtm-adam-pilot-project) or constructed ad-hoc by the `{admiral}` team. This is a prerequisite package for `{admiral}`.
* [{pharmaverseadam}](https://pharmaverse.github.io/pharmaverseadam/)---this contains test ADaM data automatically generated by running the ADaM `{admiral}` and TA package extensions templates on the [{pharmaversesdtm}](https://pharmaverse.github.io/pharmaversesdtm/) data.

Both these packages are developed by the `{admiral}` team, but can used across the pharmaverse as common, open-source test SDTM or ADaM data.

The following packages are also useful when working with ADaM datasets:

* [{metacore}](https://atorus-research.github.io/metacore/) and [{metatools}](https://pharmaverse.github.io/metatools/)---these enable users to manipulate and work with dataset metadata.
* [{xportr}](https://atorus-research.github.io/xportr/)---this provides functionality to get xpt files ready for transport.

## `{admiral}` Manifesto

For `{admiral}` and all extension packages, we prioritize providing our users with a __simple to adopt__ toolkit
that enables them to produce __readable__ and __easily constructible__ ADaM programs. The following explains
our philosophy, which we try to adhere to across the `{admiral}` family of packages.
There isn't always a clear single, straightforward rule, but there are guiding principles we adhere to for `{admiral}`.
This manifesto helps show the considerations of our developers when making decisions.

We have four design principles to achieve the main goal:

### Usability

All `{admiral}` functions should be easy to use.

* Documentation is an absolute priority. Each function reference page should cover the purpose, descriptions of each argument with permitted values, the expected input and output, with clear real-life examples---so that users don't need to dig through code to find answers.
* Vignettes that complement the functional documentation to help users see how best the functions can be applied to achieve ADaM requirements.
* Functions should be written and structured in a way that users are able to read, re-use or extend them for study specific purposes if needed (see Readability below).

### Simplicity

All `{admiral}` functions have a clear purpose.

* We try not to ever design single functions that could achieve numerous very different derivations. For example if you as a user pick up a function with >10 different arguments then chances are it is going to be difficult to understand if this function could be applied for your specific need. The intention is that arguments/parameters can influence how the output of a function is calculated, but not change the purpose of the function.

* We try to combine similar tasks and algorithms into one function where applicable to reduce the amount of repetitive functions with similar algorithms and to group together similar functionality to increase usability (e.g. one study day calculation rather than a function per variable).

* We strive to design functions that are not too general and trying to fulfill multiple, complex purposes.

* Functions should not allow expressions as arguments that are used as code snippets in function calls.

* We recommend to avoid copy and paste of complex computational algorithms or repetitive code like checks and advise to wrap them into a function. However we would also like to avoid multi-layered functional nesting, so this needs to be considered carefully to keep the nesting of 3-4 functions an exception rather than the rule.

### Findability

All `{admiral}` functions are easily findable.

* In a growing code base, across a family of packages, we make every effort to make our functions easily findable.
* We use consistent naming conventions across all our functions, and provide vignettes and ADaM templates that help users to get started and build familiarity. Each `{admiral}` family package website is searchable.
* We avoid repetitive functions that will do similar tasks (as explained above with study day example).
* Each package extension is kept focused on the specific scope, e.g. features that are relevant across multiple extension packages will be moved to the core `{admiral}` package.

### Readability

All `{admiral}` functions follow the [Programming Strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html)
that all our developers and contributors must follow, so that all our code has a high degree of consistency and readability.

* We encourage use of tidyverse (e.g. dplyr) over similar functionality existing in base R.
* For sections of code that perform the actual derivations (e.g. besides assertions or basic utilities), we try to limit nesting of too many dependencies or functions.
* Modularity is a focus---we don't try to achieve too many steps in one.
* All code has to be well commented.
* We recognize that a user or a Health Authority reviewer may have the wish to delve into the code base (especially given this open source setting), or users may need to extend/adapt the code for their study specific needs. We therefore want any module to be understandable to all, not only the `{admiral}` developers.

## References and Documentation

* Please go to [Get Started](https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html) section to start using `{admiral}`.
* Please see the [pharmaverse YouTube channel](https://www.youtube.com/channel/UCxQFEv8HNqM01DXzdQLCy6Q) for videos related to `{admiral}`.
* Please see the [Programming Strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html) to understand how functions are created.
* Please see the [FAQ](https://pharmaverse.github.io/admiral/cran-release/articles/faq.html) for the most frequent questions.
* Please see the [Contribution Model](https://pharmaverse.github.io/admiral/cran-release/CONTRIBUTING.html) for how to get involved with making contributions.
* Please see [FAQ: R and Package Versions](https://pharmaverse.github.io/admiral/cran-release/articles/faq.html#why-do-we-use-a-certain-r-version-and-package-versions-for-development) for why we develop with certain R and package versions.

## Pharmaverse Blog

If you are interested in R and Clinical Reporting, then visit the [pharmaverse blog](https://pharmaverse.github.io/blog/). This contains regular, bite-sized posts showcasing how `{admiral}` and other packages in the pharmaverse can be used to realize the vision of full end-to-end Clinical Reporting in R.

We are also always looking for keen `{admiral}` users to publish their own blog posts about how they use the package. If this could be you, feel free make an issue in the [GitHub repo](https://github.com/pharmaverse/blog) and get started!

## Recent Conference Presentations

For a full collection of `{admiral}` conference presentations over the years, please travel to our [Presentation Archive](https://pharmaverse.github.io/admiraldiscovery/articles/presentation_archive.html).

## Contact

We use the following for support and communications between user and developer community:

* [Slack](https://pharmaverse.slack.com/)---for informal discussions, Q\&A and building our user community. If you don't have access, use this [link](https://join.slack.com/t/pharmaverse/shared_invite/zt-yv5atkr4-Np2ytJ6W_QKz_4Olo7Jo9A) to join the pharmaverse Slack workspace.
* [GitHub Issues](https://github.com/pharmaverse/admiral/issues)---for direct feedback, enhancement requests or raising bugs.

## Acknowledgments

Along with the authors and contributors, thanks to the following people for their work on the package:

Jaxon Abercrombie, Mahdi About, Teckla Akinyi, Anthony Arroyo, Alex Assuied, James Black, Claudia Carlucci, Asha Chakma, 
Liming Clark, Bill Denney, Kamila Duniec, Alice Ehmann, Romain Francois, G Gayatri, Ania Golab, Alana Harris, Declan Hodges,
Solveig Holmgaard, Anthony Howard, Shimeng Huang, Samia Kabi, Leena Khatri, James Kim, John Kirkpatrick, Robin Koeger, Konstantina Koukourikou, 
Dinakar Kulkarni, Pavan Kumar, Pooja Kumari, Shan Lee, Wenyi Liu, Sadchla Mascary, Iain McCay, Jack McGavigan, Jordanna Morrish, Syed Mubasheer,
Kirill Muller, Thomas Neitmann, Yohann Omnes, Barbara O'Reilly, Celine Piraux, Hamza Rahal, Nick Ramirez, Tom Ratford, 
Sukalpo Saha, Tamara Senior, Sophie Shapcott, Vladyslav Shuliar, Eric Simms, Daniel Sjoberg, Ondrej Slama, Andrew Smith, Daniil Stefonishin,
Vignesh Thanikachalam, Michael Thorpe, Steven Ting, Ojesh Upadhyay, Franciszek Walkowiak, Enki Wang, Phillip Webster, 
Annie Yang, Andrii Yurovskyi, Junze Zhang, Kangjie Zhang, Zelos Zhu
