#  Create data/admiral_adlb.rda

# This is a TWO-step process.

library(admiral)

# First, generate script from template
adam_name="adlb"
save_path =  paste0("./data-raw/admiral_", adam_name, ".R")

use_ad_template(adam_name = adam_name,
                save_path =  save_path,
                open = F,
                overwrite=T)


# Second, run the script (manually)

