# FAQ

##### **What is `{admiral}`**?

- Think of [admiral](https://pharmaverse.github.io/admiral/) as a
  toolbox of modular blocks (R functions) to create analysis
  derivations:
  - each block has a stand alone purpose (each function provides a
    specific functionality)
  - Data Scientists can create their own blocks (create own R
    functions).
- Constructing ADaM dataset should become like building out of blocks
  that are based on [admiral](https://pharmaverse.github.io/admiral/)
  modular functions and user created modular functions.

##### **Why is the `{admiral}` logo a red butterfly?**?

- The logo is inspired from the [Red
  Admiral](https://en.wikipedia.org/wiki/Vanessa_atalanta), also known
  as *Vanessa atalanta*, which is a species of butterfly. However, our
  [extension packages](https://pharmaverse.org/e2eclinical/adam/) have
  branched out from this and have differently colored butterflies in
  their logos.

##### **Why** did we decide to start **`{admiral}`**?

- Data analysis challenges in clinical trials vary depending on
  scientific goals, therapeutic areas, indications, data sources and
  data quality. We all face the same challenge so why limit ourselves
  only to company-level adoption and crowd-sourcing to create ADaM
  datasets?
- In an effort to build ADaMs via collaboration and co-creation.
- We believed early engagement with other like-minded companies moving
  towards R could lead to our solution being shared open source as a
  framework for contribution across-industry.
- By building ADaMs like a modular building blocks, everyone can
  contribute and each module has a clear input and output to enable
  re-usable solutions.
- Users can “slot in” their own modules to address specific
  company/TA/Molecule/Study requirements.
- TA specific requirements can be open sourced again and transformed
  into a common ADaM approach for such analysis.
- We saw the long-term potential of a consistent way of producing ADaM
  and a wider community of across-industry developers contributing to
  grow the codebase to cover the infinite array of possibilities.
- We wanted to offer contributors an option to make a name for
  themselves in the Pharma open-source community and an avenue to
  collaborate with other like-minded people across the industry.
- Finally, imagine if ADaMs are built in a consistent manner with the
  same code from openly maintained functions and its impact on the
  Health Authorities, readable code, QC, talent flow…!

##### Why did we use **R as a programming language**?

- R is not an isolated software product, everyone can contribute (open
  source principle).
- People from University/Statistical talent pipeline are more likely to
  come through with R skills rather than a proprietary language .
- There is a strong data science/analytics R community.
- FDA is open to accepting R submissions and are heavy users themselves.
- Top of the line visualization/graphics - R-Shiny for interactive data
  displays and also R Markdown offers great report writing
  functionality.
- R is very popular among statisticians so new statistical methods are
  likely implemented in R before any other language.
- There might be equally suited programming languages out there -
  however at some stage we had to make a decision! :)

##### `{admiral}` offers a **toolbox of functions to facilitate ADaM**. What does that mean?

- Functions are usually parameter driven:
  - e.g. the
    [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/v1.5.0/reference/derive_vars_aage.md)
    has a parameterized start and end-date and a unit.
  - Depending on the parameters results may vary as does the
    specification.
  - Functions serve as a toolbox so the user can create their ADaM
    according to the requirements.
  - The principles, programming strategy and documentation of
    [admiral](https://pharmaverse.github.io/admiral/) are considered as
    a framework for users to contribute.

##### How does a user **know what a function does** exactly?

- Function details and its purpose, the requirements, parameters,
  dependencies and examples are documented in the header of each
  function.
- Complex functions potentially have a vignette on the
  [admiral](https://pharmaverse.github.io/admiral/) homepage to provide
  more details.
- [admiral](https://pharmaverse.github.io/admiral/) does not provide a
  link to an explicit specification in the define.xml.

##### Where can a user quickly find some **references or advice** to use a function?

- For more detailed description, please refer to the [Reference
  section](https://pharmaverse.github.io/admiral/cran-release/reference/index.html).
- [admiraldiscovery](https://pharmaverse.github.io/admiraldiscovery/articles/reactable.html)
  provides a way to look up common ADaM variables with recommended
  [admiral](https://pharmaverse.github.io/admiral/) functions to
  complete the derivation.
- Check out the [admiral Cheat
  Sheet](https://github.com/pharmaverse/admiral/blob/main/inst/cheatsheet/admiral_cheatsheet.pdf)
  as well!

##### Would `{admiral}` **create a whole** ADaM dataset?

- [admiral](https://pharmaverse.github.io/admiral/) is meant as a
  toolbox to enable Data Scientists to build ADaMs according to their
  varying analysis needs.
- [admiral](https://pharmaverse.github.io/admiral/) is not meant as a
  “click a button, out comes your ADaM” tool.
- On the [admiral](https://pharmaverse.github.io/admiral/) website,
  example scripts are provided which can be used as a starting point to
  create an ADaM (see at the end of a vignette).

##### In **which order** does a user need to execute the functions?

- Guidance is available for ADSL, BDS and OCCDS ADaM structure including
  template scripts.

##### Is the `{admiral}` package **validated**?

- All functions are reviewed and tested (see \[What is provided around
  **function testing**?\]) to ensure that they work as described in the
  documentation.
- Test cases for each function are be part of the R package.
- Users can add to the tests or provide additional feedback.
- The testing the [admiral](https://pharmaverse.github.io/admiral/) team
  does for each function does **not replace the QC and validation
  process at each company**.
- A [GitHub
  action](https://github.com/marketplace/actions/r-package-validation-report)
  (using open source packages) exists to generate a validation report
  for an R package, which would be an option for any company to use. An
  example report using an earlier version of
  [admiral](https://pharmaverse.github.io/admiral/) exists.
  [here](https://github.com/insightsengineering/thevalidatoR/blob/main/readme_files/report-0.1-admiral.pdf)
  as an illustration.

##### What is available around **function testing**?

- Unit tests for reliability of each function - available as part of
  open source release.
- Some integration testing is be done to ensure functions can be called
  together to create ADaM (e.g. even via the internal testing teams).
- Guidance for testing and documentation expectations of community
  contribution functions. Then it is for each company to cover the
  following:
  - validation to be able to use the package on company-specific SCE for
    GxP purposes and associated audit evidence
  - strategy of how the use of
    [admiral](https://pharmaverse.github.io/admiral/) fits into
    company-specific quality assurance process (double programming
    comparison versus your company-specific legacy ADaM solution could
    be appropriate until confidence builds)
  - see our guidance on [unit
    testing](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html).

##### Does **admiral provide harmonized define.xml** or submittable specifications for functions?

- No. The functions are documented via programming headers, the
  define.xml is the responsibility of the end user.
- Functions are usually generalized and not specific. (see [Admiral
  offers a **toolbox of functions to facilitate ADaM**. What does that
  mean?](#admiral-offers-a-toolbox-of-functions-to-facilitate-adam--what-does-that-mean))
- The users are responsible to make sure they use the functions and
  their parameters in the right way to ensure alignment with their
  define.xml

##### Does `{admiral}` provide ADaM IG **CDISC compliant** datasets?

- Although [admiral](https://pharmaverse.github.io/admiral/) follows
  CDISC standards it does not claim that the dataset resulting from
  calling [admiral](https://pharmaverse.github.io/admiral/) functions is
  ADaM compliant. This has to be ensured by the user.

##### How much of the **ADaM IG is covered by `{admiral}`?**

- ADaM IG is a standard framework without a specific number of datasets
  or variables, so it cannot be used as a specific baseline to answer
  that question.
- We provide guidance for each ADaM dataset structure (ADSL, OCCDS and
  BDS) that highlights which functionality
  [admiral](https://pharmaverse.github.io/admiral/) covers. (see [In
  **which order** does a user need to execute the
  functions?](#in-which-order-does-a-user-need-to-execute-the-functions)).
- The guidance also highlights the gaps to be filled by the user
  (e.g. timing, ranges).

##### Is there be a user/**contribution** guide?

- Our [programming
  strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html)
  serves as a framework for users how to create their own functions.
- Please see the [contribution
  model](https://pharmaverse.github.io/admiral/cran-release/CONTRIBUTING.html)
  for a detailed description on how to contribute to the project.

##### Do `{admiral}` functions **preserve variable attributes** like labels?

- Most [admiral](https://pharmaverse.github.io/admiral/) functions don’t
  preserve variable attributes because they rely on
  [dplyr](https://dplyr.tidyverse.org) operations, which generally don’t
  preserve attributes.
- For [admiral](https://pharmaverse.github.io/admiral/) functions which
  preserve attributes like
  [`convert_blanks_to_na()`](https:/pharmaverse.github.io/admiral/v1.5.0/reference/convert_blanks_to_na.md)
  and
  [`convert_na_to_blanks()`](https:/pharmaverse.github.io/admiral/v1.5.0/reference/convert_na_to_blanks.md)
  this is explicitly mentioned in the documentation.
- The recommended approach is to apply variable labels and other
  metadata as a final step in your data derivation process using
  packages like [metacore](https://atorus-research.github.io/metacore/),
  [metatools](https://pharmaverse.github.io/metatools/), and
  [xportr](https://atorus-research.github.io/xportr/).

##### Are there any **presentations** available about `{admiral}`?

- For a full collection of
  [admiral](https://pharmaverse.github.io/admiral/) conference
  presentations over the years, please travel to our [Presentation
  Archive](https://pharmaverse.github.io/admiraldiscovery/articles/presentation_archive.html).
