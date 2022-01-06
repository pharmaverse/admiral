
library(tibble)
library(stringr)
system('curl -H "Accept: application/vnd.github.v3+json" https://api.github.com/repos/Roche-GSK/admiral/issues > issues.txt')

issues <- readLines("issues.txt")

issues_df <- enframe(unlist(issues))

issues_filter <- issues_df %>% filter(str_detect(value, "number"))
