# Contribution to {admiral} 

This outlines how to propose a change to the admiral package. For more detailed info about contributing to {admiral}, and other [pharmaverse packages](https://pharmaverse.org/), please see the [Contribution Guide](https://pharmaverse.github.io/admiral/cran-release/CONTRIBUTING.html) as well as other Developer Guides in the Articles section of the [{admiraldev} website](https://pharmaverse.github.io/admiraldev/).

Please note that we try to align to best practices used in other R packages' development processes - so veteran developers should be familiar with our processes. However, we do deviate slightly from some best practices and we advise all new contributors to review our package documentation accordingly.


# Basics of Contribution
 🦋 For each new contribution, the user creates an issue on the issue tab on [GitHub](https://github.com/pharmaverse/admiral/issues) to put it in our backlog. The issues can range from bug identification and/or fixes, enhancements to functions, documentation, tests or new features.   

 🦋 We advise you to contact us when an [issue](https://github.com/pharmaverse/admiral/issues) is created via [Slack](https://app.slack.com/client/T028PB489D3/C02M8KN8269) (If you don't have access, use this [link](https://join.slack.com/t/pharmaverse/shared_invite/zt-yv5atkr4-Np2ytJ6W_QKz_4Olo7Jo9A) to join).  We can discuss details or align expectations if you are not familiar with the `{admiral}` philosophy and programming strategy. The team will try to review the issues within the next backlog meeting and give some initial feedback. Since we are not a 100% fully resourced software development team it might be that some issues will take longer to respond to depending on the amount of overall issues. 

 🦋 We advise you to familiarize yourself with our [programming strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html), guidance for [GitHub usage](https://pharmaverse.github.io/admiraldev/articles/git_usage.html) and [unit testing](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html).

 🦋 All newly [created issues](https://github.com/pharmaverse/admiral/issues) will be reviewed within the next backlog meeting and the creator will receive an initial feedback via a comment. Someone from the core development team will then triage new issues by assigning the appropriate labels (such as “user request” so we can easily identify new requests).
 
 🦋 Issues are meant to be taken on by users from the Pharma programming community and not exclusively by the `{admiral}` team from GSK and Roche.

# Contribution Model

## Type 1: Contribution without Code 

  🦋  User creates an issue and ideally contacts an `{admiral}` team member via [Slack](https://app.slack.com/client/T028PB489D3/C02M8KN8269).
  
  🦋  The `{admiral}` core development team will contact the issue creator as soon as possible to discuss further details.
 

## Type 2: Contribution with Code

  🦋  First, the user creates an issue or comments on an existing issue to notify that they’d like to contribute code.
  
  🦋  Follow our development process step-by-step guide.
  
  🦋  We advise to contact an `{admiral}` core development team directly via [Slack](https://app.slack.com/client/T028PB489D3/C02M8KN8269) before submitting code for complex functionality.

## Detailed Development Process

If you decide to contribute with code and you're ready to make your first code contribution, this detailed development process step-by-step guide will help tie all the other detailed vignettes together to give you the simplest experience of helping to grow and enhance our codebase.

1.  Create a new feature branch from the `main` branch
    following the naming convention and pull the latest changes - as
    detailed on the [GitHub
    usage](https://pharmaverse.github.io/admiraldev/articles/git_usage.html#working-with-feature-branches-1) guide.
2.  Familiarize yourself with the `{admiral}` [programming
    strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html), and then make the required
    code updates.
3.  Before making a pull request, check the [Pull Request Review
    Guidance](https://pharmaverse.github.io/admiraldev/articles/pr_review_guidance.html) & the following checklist of
    common things developers miss:
    a.  Is all your code formatted according to the
        [tidyverse](https://style.tidyverse.org/) style guide?
    b.  Did you create/add appropriate [unit
        tests](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html#writing-unit-tests-in-admiral)?
    c.  If you removed/replaced any function and/or function parameters,
        did you fully follow the [deprecation
        guidance](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#deprecation)?
    d.  Did you update the
        [documentation]https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#function-header-documentation)?
        If so, remember to run `devtools::document()` and include the
        updated `NAMESPACE` and `.Rd` files in `man/`.
    e.  Does your code update have any impact on the [ADaM
        template](https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#starting-a-script)
        R scripts stored in `inst/templates`?
    f.  Does your code update have any impact on the vignettes stored in
        vignettes?
    g.  Did you update the Changelog `NEWS.md`?
    h.  Did you build `{admiral}` site `pkgdown::build_site()` and check
        that all affected examples are displayed correctly and that all
        new functions occur on the
        "[Reference](https://pharmaverse.github.io/admiral/cran-release/reference/)" page?
4.  Once happy with all the updates, make a [pull
    request](https://pharmaverse.github.io/admiraldev/articles/git_usage.html#pull-request) to merge to the `main` branch
    and link the issue so that it closes after successful merging.
5.  Check that there are no merge conflicts. If there are any, fix them
    before requesting review. See [solving merge
    conflicts](https://pharmaverse.github.io/admiraldev/articles/git_usage.html#solving-merge-conflicts-in-the-terminal-on-rstudio)
    guidance.
6.  Check the results of the automated `R-CMD check` and `lintr` checks
    and if any issues consult this
    [guide](https://pharmaverse.github.io/admiraldev/articles/pr_review_guidance.html#common-r-cmd-check-issues).
7.  Assign a reviewer from the `{admiral}` core development team - this
    could be anyone you discussed the issue with previously via Slack or
    GitHub. If unsure, add a comment that the pull request is ready for
    review and add the `@pharmaverse/admiral` tag to it.
8.  Once the review is completed, the reviewer will merge the PR and
    this will then automatically delete the feature branch.

Finally, just a note to say from the core developers that we hugely
appreciate you taking the time to contribute to `{admiral}`. Don't be
offended if during review we send requests back to you, as the
expectations are high so that we can ensure the `{admiral}` codebase
remains robust and consistent. The best way to learn here is always to
jump in and get involved, so please don't be afraid you'll make mistakes
along the way -- we all have and continue to do so, and that's what the
reviews are for. Also if ever you get stuck don't hesitate to reach out
for support via the [Slack
channel](https://pharmaverse.slack.com/).
 ***Welcome to our `{admiral}` community!***
