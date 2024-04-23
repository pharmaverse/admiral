# list_all_templates Test 2: Error Message is returned if package is not installed

    Code
      list_all_templates(package = "non-existing-package")
    Condition
      Error in `list_all_templates()`:
      ! No package called non-existing-package is installed and hence no templates are available.

# use_ad_template Test 4: Error Message is returned if no ADaM template is available

    Code
      suppressMessages(use_ad_template("adxx", save_path = file, open = FALSE))
    Condition
      Error in `use_ad_template()`:
      ! No template for ADXX available in package admiral.
      i Run `admiral::list_all_templates("admiral")` to get a list of all available ADaM templates.

